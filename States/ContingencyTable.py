"""
    FeatureCloud SPLINK Application
    Copyright 2021 Mohammad Bakhtiari. All Rights Reserved.
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
        http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
"""

from FeatureCloud.app.engine.app import app_state, AppState, Role, SMPCOperation
import numpy as np
import multiprocessing
import threading
from Utils.gwas_dataset import PhenotypeValue, SnpValue
from SplinkStates.client import SplinkClient
from SplinkStates.server import SplinkServer
from Utils.utils import share_attrs, load_attrs
from States import ALGORITHM


@app_state('Contingency_Table', Role.BOTH)
class ContingencyTable(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_Contingency_Table', Role.COORDINATOR)
        self.register_transition('terminal', Role.PARTICIPANT)

    def run(self) -> str or None:
        load_attrs(self)
        data = self.await_data()
        global_minor_allele_names, global_major_allele_names, data = data
        self.minor_allele(global_minor_allele_names, global_major_allele_names)
        if self.load('config')['algorithm'] in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]:
            self.snp_indices = data
            contingency_tables = self.contingency_table_step()
            self.send_data_to_coordinator(data=contingency_tables, use_smpc=self.load('smpc_used'))
        elif self.load('config')['algorithm'] == ALGORITHM.LOGISTIC_REGRESSION:
            pass
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_Contingency_Table'
        return 'terminal'

    def minor_allele(self, global_minor_allele_names, global_major_allele_names):

        for snp_index in global_minor_allele_names.keys():
            # if local minor/major allele is different from the global minor/major allele
            if self.second_allele_names[snp_index] != global_major_allele_names[snp_index]:
                # swap the local minor and major allele names
                self.first_allele_names[snp_index] = global_minor_allele_names[snp_index]
                self.second_allele_names[snp_index] = global_major_allele_names[snp_index]

                # inverse the mapping of the SNP values 0 -> 2 and 2 -> 0
                self.snp_values[snp_index] = np.where(self.snp_values[snp_index] == 2, -3, self.snp_values[snp_index])
                self.snp_values[snp_index] = np.where(self.snp_values[snp_index] == 0, 2, self.snp_values[snp_index])
                self.snp_values[snp_index] = np.where(self.snp_values[snp_index] == -3, 0, self.snp_values[snp_index])
        self.attrs_to_share += ['snp_values', 'first_allele_names', 'second_allele_names']

    # ##### contingency table step related functions
    def contingency_table_step(self):
        """ Compute local contingency table for the chunk"""
        # get SNP indices for which contingency table should be computed

        self.current_chunk_size = len(self.snp_indices)

        # queue
        queue_contingency_tables = multiprocessing.Queue()

        # thread to read from the queues
        thread_read_contingency_tables_queue = threading.Thread(target=self.read_queue_contingency_tables,
                                                                args=(queue_contingency_tables,))
        thread_read_contingency_tables_queue.daemon = True
        thread_read_contingency_tables_queue.start()

        # processes to compute the local contingency tables for the sub-chunks
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices, self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_contingency_table,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    queue_contingency_tables,))

            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        thread_read_contingency_tables_queue.join()

        # close queue
        queue_contingency_tables.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionaries to lists;
        contingency_tables = list()
        for snp_index in sorted(self.snp_indices):
            contingency_tables.append(self.contingency_tables[snp_index])
        return contingency_tables

    def read_queue_contingency_tables(self, queue_contingency_tables):
        while len(self.contingency_tables) < self.current_chunk_size:
            cont_table = queue_contingency_tables.get()
            self.contingency_tables.update(cont_table)

        # compute contingency table for a set of SNPs

    def compute_contingency_table(self, start_snp_index, end_snp_index, contingency_table_queue):
        """ Compute local contingency table for a sub-chunk """

        contingency_tables = dict()

        for snp_index in np.arange(start_snp_index, end_snp_index):
            if snp_index not in self.snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                contingency_table_queue.put(contingency_tables)
                contingency_tables = dict()

            minor_allele = self.first_allele_names[snp_index]
            major_allele = self.second_allele_names[snp_index]

            # minor-case
            minor_case_count = self.compute_allele_count(snp_index, minor_allele, PhenotypeValue.CASE)

            # major-case
            major_case_count = self.compute_allele_count(snp_index, major_allele, PhenotypeValue.CASE)

            # minor-control
            minor_control_count = self.compute_allele_count(snp_index, minor_allele, PhenotypeValue.CONTROL)

            # major-control
            major_control_count = self.compute_allele_count(snp_index, major_allele, PhenotypeValue.CONTROL)

            # contingency table value: [minor-case, major-case, minor-control, major-control]
            contingency_tables[snp_index] = np.array([minor_case_count,
                                                      major_case_count,
                                                      minor_control_count,
                                                      major_control_count])

        # put the remaining contingency tables into queue
        contingency_table_queue.put(contingency_tables)

    def compute_allele_count(self, snp_index, allele_name, trait):
        """ Compute allele count for minor-case, minor-control, major-case, and major-control """
        x_matrix, phenotype_values = self.get_x_matrix_y_vector(snp_index,
                                                                covariate=self.load('config')['covariates'],
                                                                algorithm=self.load('config')['algorithm'])
        snp_values = x_matrix[:, 1]

        trait_indices = np.where(phenotype_values == trait)[0]
        trait_snp_values = snp_values[trait_indices]

        if allele_name == self.first_allele_names[snp_index]:
            return int(2 * np.where(trait_snp_values == SnpValue.HOMOZYGOTE_00)[0].size +
                       np.where(trait_snp_values == SnpValue.HETEROZYGOTE)[0].size)

        if allele_name == self.second_allele_names[snp_index]:
            return int(2 * np.where(trait_snp_values == SnpValue.HOMOZYGOTE_11)[0].size +
                       np.where(trait_snp_values == SnpValue.HETEROZYGOTE)[0].size)


@app_state('Aggregate_Contingency_Table', Role.COORDINATOR)
class AggregateContingencyTable(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Non_Missing_counts', Role.COORDINATOR)
        self.register_transition('terminal', Role.COORDINATOR)

    def run(self) -> str or None:
        load_attrs(self)
        contingency_tables = self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
        self.contingency_table_step(contingency_tables)
        # if this is not the last chunk, set up the next chunk of SNPs
        if not self.is_last_chunk():
            self.load_attrs()
            data_to_send = self.setup_next_chunk()
            self.broadcast_data(data_to_send)
            share_attrs(self)
            return 'Non_Missing_counts'
        # if this is the last chunk, generate the manhattan plot first, and then, tell clients to download the results
        self.manhattan_plot(manhatan_filename=self.load('output_files')['manhattan'][0], log=self.log)
        # self.set_step(HyFedProjectStep.RESULT)
        share_attrs(self)
        return 'terminal'

    def contingency_table_step(self, contingency_tables):
        """ Compute global chi-square, odd ratio, and p-values using the aggregated contingency tables """

        # convert global contingency table from list to dictionary
        snp_counter = -1
        for snp_index in sorted(self.considered_snp_indices):
            snp_counter += 1
            self.contingency_tables[snp_index] = contingency_tables[snp_counter]
        self.attrs_to_share += ['contingency_tables']
        # compute the results (i.e. MAF, chi-square, odd-ratio, and p-values) for the chunk
        self.compute_results_chi_square()

        # add chromosome number, base pair distance, and p-value of the current chunk to results for all chunks
        self.append_to_results_all_chunks()
        self.log(self.load('output_files')['output'][0])
        # save the results using a separate process
        save_process = multiprocessing.Process(target=self.save_results_chi_square,
                                               args=(self.load('output_files')['output'][0],))
        save_process.daemon = True
        save_process.start()
        save_process.join()
        save_process.terminate()

    def compute_results_chi_square(self):
        """ Compute MAF for case/control, chi-square, odd-ratio, and p-values for chi-square algorithm """
        msg = self.compute_maf()
        self.log(msg)
        msg = self.compute_chi_square_values()
        self.log(msg)
        msg = self.compute_odd_ratio_values()
        self.log(msg)
        msg = self.compute_p_values(self.load('config')['algorithm'], self.load('config')['covariates'])
        self.log(msg)
