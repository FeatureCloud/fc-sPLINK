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

from FeatureCloud.app.engine.app import app_state, AppState, Role, SMPCOperation, LogLevel
from FeatureCloud.app.engine.app import State as op_state
from time import sleep
import numpy as np
import multiprocessing
import threading
from Utils.gwas_dataset import SnpValue
from SplinkStates.client import SplinkClient
from SplinkStates.server import SplinkServer
from Utils.utils import share_attrs, load_attrs
from States import ALGORITHM


@app_state('Non_Missing_counts', Role.BOTH)
class NonMissingCounts(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_Non_Missing_Counts', Role.COORDINATOR)
        self.register_transition('Contingency_Table', Role.PARTICIPANT)
        self.register_transition('Beta_Linear', Role.PARTICIPANT)
        self.register_transition('Beta_Logistic', Role.PARTICIPANT)

    def run(self) -> str or None:
        """
        Sharing:
        non_missing_sample_counts
        allele_counts

        Returns
        -------

        """
        load_attrs(self)
        self.current_chunk, self.total_chunks, self.snp_indices, self.chunk_start_index, self.chunk_end_index = \
            self.await_data()
        self.current_chunk_size = len(self.snp_indices)
        non_missing_sample_counts, allele_counts = self.non_missing_counts()
        self.send_data_to_coordinator(data=non_missing_sample_counts, use_smpc=self.load('smpc_used'))
        if self.is_coordinator:
            # self.non_missing_sample_counts = \
            #     self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
            self.store('non_missing_sample_counts',
                       self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
                       )
        else:
            sleep(5)
        self.send_data_to_coordinator(data=allele_counts, use_smpc=self.load('smpc_used'))


        # if self.load('config')['algorithm'] in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]:
        #     self.send_data_to_coordinator(data=allele_counts, use_smpc=self.load('smpc_used'))
        # else:
        #     beta_values = self.get_beta_values()
        #     # share initial beta values (excluding ignored SNPs) for which gradient, Hessian, and log likelihood
        #     # should be computed in clients
        #     # self.send_data_to_coordinator(data=[beta_values, self.current_beta_iteration])
        #     self.send_data_to_coordinator(data=beta_values)
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_Non_Missing_Counts'
        if self.load('config')['algorithm'] == ALGORITHM.CHI_SQUARE:
            return 'Contingency_Table'
        if self.load('config')['algorithm'] == ALGORITHM.LINEAR_REGRESSION:
            return 'Beta_Linear'
        return 'Beta_Logistic'

    def non_missing_counts(self):
        self.set_sub_chunk_indices(self.chunk_start_index, self.chunk_end_index)
        # init count dictionaries
        self.non_missing_sample_counts, self.allele_counts = {}, {}

        # queues
        queue_non_missing = multiprocessing.Queue()
        queue_allele_counts = multiprocessing.Queue()

        # threads to read from the queues
        thread_read_non_missing_queue = threading.Thread(target=self.read_queue_non_missing, args=(queue_non_missing,))
        thread_read_non_missing_queue.daemon = True
        thread_read_non_missing_queue.start()

        thread_read_allele_counts_queue = threading.Thread(target=self.read_queue_allele_counts,
                                                           args=(queue_allele_counts,))
        thread_read_allele_counts_queue.daemon = True
        thread_read_allele_counts_queue.start()

        # start processes to compute the local non-missing sample counts as well as first/second allele counts for sub-chunks
        process_list = []
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices, self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_non_missing_counts,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    queue_non_missing, queue_allele_counts,))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read threads to be done
        thread_read_non_missing_queue.join()
        thread_read_allele_counts_queue.join()

        # close queues
        queue_non_missing.close()
        queue_allele_counts.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionaries to lists;
        # IMPORTANT: sorted(self.snp_indices) should always be used to ensure the order between list and set related SNP indices

        # non_missing_sample_counts, allele_counts = [], []
        # for snp_index in sorted(self.snp_indices):
        #     non_missing_sample_counts.append(self.non_missing_sample_counts[snp_index])
        #     allele_counts.append(self.allele_counts[snp_index])
        non_missing_sample_counts = [self.non_missing_sample_counts[ind] for ind in self.snp_indices]
        allele_counts = [self.allele_counts[ind] for ind in self.snp_indices]
        # python list of scalars must be converted to a numpy array if compensator flag is set
        non_missing_sample_counts = np.array(non_missing_sample_counts)
        return non_missing_sample_counts, allele_counts

    # ##### multi-processing functions
    def set_sub_chunk_indices(self, start_snp_index, end_snp_index):
        """ Determine start/end indices for sub-chunks assigned to each process/core """

        if end_snp_index <= start_snp_index:
            self.log("end_snp_index must be greater than start_snp_index!", LogLevel.FATAL)
            self.update(state=op_state.ERROR)

        # ensure each process/core will compute at least one SNP statistics
        if self.current_chunk_size < self.load('config')['cpu_cores']:
            cpu_cores = 1
        else:
            cpu_cores = self.load('config')['cpu_cores']

        sub_chunk_size = int(np.ceil(self.current_chunk_size / cpu_cores))

        start_indices = np.arange(start_snp_index, end_snp_index, sub_chunk_size)
        end_indices = start_indices + sub_chunk_size
        end_indices[-1] = end_snp_index
        self.sub_chunk_start_indices = start_indices
        self.sub_chunk_end_indices = end_indices

    def compute_non_missing_counts(self, start_index, end_index, queue_non_missing, queue_allele_counts):
        """ Compute local non-missing sample count as well as first/second allele count for a sub-chunk """

        # init dictionaries
        non_missing_sample_counts, allele_counts = {}, {}
        for snp_index in np.arange(start_index, end_index):

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_non_missing.put(non_missing_sample_counts)
                queue_allele_counts.put(allele_counts)
                non_missing_sample_counts, allele_counts = {}, {}

            # non-missing sample count
            x_matrix, y_vector = self.get_x_matrix_y_vector(snp_index,
                                                            covariate=self.load('config')['covariates'],
                                                            algorithm=self.load('config')['algorithm'])
            non_missing_sample_counts[snp_index] = y_vector.size

            # allele count
            # snp_values = self.load('snp_values')[snp_index]
            snp_values = self.snp_values[snp_index]
            first_allele_count = int(
                2 * np.where(np.array(snp_values) == SnpValue.HOMOZYGOTE_00)[0].size +
                np.where(np.array(snp_values) == SnpValue.HETEROZYGOTE)[0].size
            )

            second_allele_count = int(
                2 * np.where(np.array(snp_values) == SnpValue.HOMOZYGOTE_11)[0].size +
                np.where(np.array(snp_values) == SnpValue.HETEROZYGOTE)[0].size
            )

            # to stick with the correct mapping allele_name -> allele count, which is based on the sorted allele names in the server
            # if self.load('first_allele_names')[snp_index] < self.load('second_allele_names')[snp_index]:
            if self.first_allele_names[snp_index] < self.second_allele_names[snp_index]:
                allele_counts[snp_index] = np.array([first_allele_count, second_allele_count])
            else:
                allele_counts[snp_index] = np.array([second_allele_count, first_allele_count])

        # put remaining results in the corresponding queues
        queue_non_missing.put(non_missing_sample_counts)
        queue_allele_counts.put(allele_counts)

    # ##### Queue functions
    def read_queue_non_missing(self, queue_non_missing):
        while len(self.non_missing_sample_counts) < self.current_chunk_size:
            sample_count_non_missing = queue_non_missing.get()
            self.non_missing_sample_counts.update(sample_count_non_missing)

    def read_queue_allele_counts(self, queue_allele_counts):
        while len(self.allele_counts) < self.current_chunk_size:
            count_alleles = queue_allele_counts.get()
            self.allele_counts.update(count_alleles)

    # def get_beta_values(self):
    #     # initialize log_likelihood and beta values
    #     self.considered_in_process_snp_indices = self.considered_snp_indices.intersection(
    #         self.in_process_snp_indices)
    #     beta_values = dict()  # beta values shared with clients
    #     for snp_index in self.considered_in_process_snp_indices:
    #         # 2 for snp and intercept columns
    #         beta_values[snp_index] = np.array([0.0 for _ in range(0, len(self.load('config')['covariates']) + 2)])
    #         # beta_values[snp_index] = np.array([0.0 for _ in range(0, len(self.covariates) + 2)])
    #         self.beta_values[snp_index] = beta_values[snp_index]
    #         self.log_likelihood_values[snp_index] = None


@app_state('Aggregate_Non_Missing_Counts', Role.COORDINATOR)
class AggregateNonMissingCounts(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Contingency_Table', Role.COORDINATOR)
        self.register_transition('Beta_Linear', Role.COORDINATOR)
        self.register_transition('Beta_Logistic', Role.COORDINATOR)


    def run(self) -> str or None:
        """
        Share:
        minor_allele_names_considered
        major_allele_names_considered
        Returns
        -------

        """
        load_attrs(self)
        allele_counts = self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
        # allele_counts = self.gather_data()
        non_missing_sample_counts = self.load('non_missing_sample_counts')
        # non_missing_sample_counts = self.non_missing_sample_counts
        data_to_send = self.aggregate_non_missing_counts(non_missing_sample_counts, allele_counts)
        # To share
        # From nonmissing count : minor_allele_names_considered, major_allele_names_considered
        #
        if self.load('config')['algorithm'] in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]:
            data_to_send.append(self.considered_snp_indices)
        else:
            beta_values, converged = self.aggregate_minor_allele()
            data_to_send.append(beta_values)
            data_to_send.append(converged)
        # self.broadcast_data(data=[minor_allele_names_considered, major_allele_names_considered, data])
        self.broadcast_data(data_to_send)
        share_attrs(self)
        if self.load('config')['algorithm'] == ALGORITHM.CHI_SQUARE:
            return 'Contingency_Table'
        if self.load('config')['algorithm'] == ALGORITHM.LINEAR_REGRESSION:
            return 'Beta_Linear'
        return 'Beta_Logistic'

    def aggregate_non_missing_counts(self, non_missing_sample_counts, allele_counts):

        # determine global minor/major allele for each SNP
        minor_allele_names = np.argmin(allele_counts, axis=1)
        major_allele_names = 1 - minor_allele_names

        # get the minor/major allele count for each SNP
        minor_allele_counts = np.min(allele_counts, axis=1)
        major_allele_counts = np.max(allele_counts, axis=1)

        # compute minor/major allele frequency for each SNP
        allele_counts_total = np.sum(allele_counts, axis=1)
        minor_allele_frequencies = minor_allele_counts / allele_counts_total
        major_allele_frequencies = major_allele_counts / allele_counts_total

        # store global non-missing sample count, minor/major allele names/counts/frequencies in the dictionaries;
        snp_counter = -1
        for snp_index in sorted(self.considered_snp_indices.copy()):
            snp_counter += 1
            self.non_missing_sample_counts[snp_index] = non_missing_sample_counts[snp_counter]
            self.allele_counts[snp_index] = allele_counts[snp_counter]
            self.minor_allele_names[snp_index] = self.allele_names[snp_index][minor_allele_names[snp_counter]]
            self.major_allele_names[snp_index] = self.allele_names[snp_index][major_allele_names[snp_counter]]
            self.minor_allele_counts[snp_index] = minor_allele_counts[snp_counter]
            self.major_allele_counts[snp_index] = major_allele_counts[snp_counter]
            self.minor_allele_frequencies[snp_index] = minor_allele_frequencies[snp_counter]
            self.major_allele_frequencies[snp_index] = major_allele_frequencies[snp_counter]
        self.attrs_to_share += ['minor_allele_names', 'major_allele_names', 'minor_allele_counts',
                                'major_allele_counts',
                                'minor_allele_frequencies', 'major_allele_frequencies']
        # share the global minor/major allele names
        minor_allele_names_considered = dict()
        major_allele_names_considered = dict()
        for snp_index in self.considered_snp_indices:
            minor_allele_names_considered[snp_index] = self.minor_allele_names[snp_index]
            major_allele_names_considered[snp_index] = self.major_allele_names[snp_index]
        return [minor_allele_names_considered, major_allele_names_considered]

    def aggregate_minor_allele(self):
        if self.load('config')['algorithm'] in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]:
            return self.considered_snp_indices
        elif self.load('config')['algorithm'] == ALGORITHM.LOGISTIC_REGRESSION:
            # initialize log_likelihood and beta values
            self.considered_in_process_snp_indices = self.considered_snp_indices.intersection(
                self.in_process_snp_indices)
            beta_values = dict()  # beta values shared with clients
            for snp_index in self.considered_in_process_snp_indices:
                # 2 for snp and intercept columns
                beta_values[snp_index] = np.array([0.0 for _ in range(0, len(self.load('config')['covariates']) + 2)])
                # beta_values[snp_index] = np.array([0.0 for _ in range(0, len(self.covariates) + 2)])
                self.beta_values[snp_index] = beta_values[snp_index]
                self.log_likelihood_values[snp_index] = None

            # share initial beta values (excluding ignored SNPs) for which gradient, Hessian, and log likelihood
            # should be computed in clients
            self.beta_values = beta_values
            return beta_values, False
            # return [beta_values, self.current_beta_iteration]
