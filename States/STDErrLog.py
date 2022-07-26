"""
    FeatureCloud SPLINK Application
    Copyright 2022 Mohammad Bakhtiari. All Rights Reserved.
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


@app_state('STD_Error_Logistic', Role.BOTH)
class STDErrorLogistic(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_STD_Error_Logistic', Role.COORDINATOR)
        self.register_transition('Limbo', Role.PARTICIPANT)

    def run(self) -> str or None:
        beta_values = self.load('beta_values')
        load_attrs(self)
        hessian_matrices = self.std_error_logistic_step(beta_values)
        self.send_data_to_coordinator(hessian_matrices, use_smpc=self.load('smpc_used'))
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_STD_Error_Logistic'
        return 'Limbo'

    # #####  std error step related functions for logistic regression algorithm
    def std_error_logistic_step(self, beta_values):
        """ Compute local hessian matrices for the chunk """

        # set Hessian dictionary to empty at the beginning of the chunk
        self.hessian_matrices = dict()

        # queue
        queue_hessian = multiprocessing.Queue()

        # thread to read from the queue
        hessian_read_thread = threading.Thread(target=self.read_queue_hessian, args=(queue_hessian,))
        hessian_read_thread.daemon = True
        hessian_read_thread.start()

        # global beta values

        self.current_chunk_size = len(beta_values)

        # processes to compute the local Hessian matrices for the sub-chunks
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices,
                                                              self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_hessian_matrices,
                                              args=(start_index_sub_chunk, end_index_sub_chunk, beta_values,
                                                    queue_hessian,))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        hessian_read_thread.join()

        # close queues
        queue_hessian.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionary to list
        hessian_matrices = list()
        for snp_index in sorted(beta_values.keys()):
            hessian_matrices.append(self.hessian_matrices[snp_index])
        return hessian_matrices
        # # share noisy local Hessian matrices with the server and noise with compensator
        # self.local_parameters[SplinkLocalParameter.HESSIAN] = hessian_matrices
        # self.set_compensator_flag({SplinkLocalParameter.HESSIAN: DataType.LIST_NUMPY_ARRAY_FLOAT})

    def compute_hessian_matrices(self, start_index, end_index, beta_values, queue_hessian):
        """ Compute local Hessian matrices for a sub-chunk """

        hessian_matrices = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in beta_values.keys():
                continue

            # put results in the queues whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_hessian.put(hessian_matrices)
                hessian_matrices = dict()

            # Hessian matrix
            x_matrix, y_vector = self.get_x_matrix_y_vector(snp_index,
                                                            covariate=self.load('config')['covariates'],
                                                            algorithm=self.load('config')['algorithm'])
            beta_vector = beta_values[snp_index].reshape(-1, 1)
            x_beta_product = np.dot(x_matrix, beta_vector)
            y_predicted = 1 / (1 + np.exp(-x_beta_product))
            hessian_matrices[snp_index] = np.dot(np.multiply(x_matrix.T, (y_predicted * (1 - y_predicted)).T),
                                                 x_matrix)

        queue_hessian.put(hessian_matrices)


@app_state('Aggregate_STD_Error_Logistic', Role.COORDINATOR)
class AggregateSTDErrorLogistic(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Non_Missing_counts', Role.COORDINATOR)
        self.register_transition('terminal', Role.COORDINATOR)

    def run(self) -> str or None:
        hessian_matrices = self.aggregate_data(use_smpc=self.load('smpc_used'))
        load_attrs(self)
        self.std_error_logistic_step(hessian_matrices)
        # if this is not the last chunk, set up the next chunk of SNPs
        if not self.is_last_chunk():
            data_to_send = self.setup_next_chunk()
            self.broadcast_data(data_to_send)
            share_attrs(self)
            return 'Non_Missing_counts'
        # if this is the last chunk, generate the manhattan plot first, and then, tell clients to download the results
        self.manhattan_plot(self.load('output_files')['manhattan'][0], self.log)
        # self.set_step(HyFedProjectStep.RESULT)
        return 'terminal'

    # ##### logistic regression std-error step related functions
    def std_error_logistic_step(self, hessian_matrices):
        """ Compute logistic regression standard error values using the aggregated Hessian matrices for the chunk """

        # convert list to dictionary
        self.hessian_matrices = dict()
        snp_counter = -1
        for snp_index in sorted(self.considered_snp_indices):
            snp_counter += 1
            self.hessian_matrices[snp_index] = hessian_matrices[snp_counter]

        # initialize std_error_values as an empty dictionary
        self.std_error_values = dict()

        # queue
        queue_std_error = multiprocessing.Queue()

        # thread to read from the queue
        std_error_read_thread = threading.Thread(target=self.read_queue_std_error, args=(queue_std_error,))
        std_error_read_thread.daemon = True
        std_error_read_thread.start()

        # processes to compute the std error values for the sub-chunks
        sub_chunk_start_indices, sub_chunk_end_indices = self.get_start_end_indices(cpu_cores=8)
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(sub_chunk_start_indices, sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.calculate_std_error_logistic_sub_chunk,
                                              args=(
                                                  start_index_sub_chunk, end_index_sub_chunk, queue_std_error,))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        std_error_read_thread.join()

        # close queues
        queue_std_error.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # update ignored index set
        for snp_index in self.considered_snp_indices:
            if self.std_error_values[snp_index][0] == "NA":
                self.considered_snp_indices.discard(snp_index)
                self.t_stat_values[snp_index] = self.std_error_values[snp_index]
                self.p_values[snp_index] = self.std_error_values[snp_index]
                continue

        # compute the results (i.e. t-stats and p-values) for the chunk
        self.compute_results_regression(self.load('config')['algorithm'], self.load('config')['covariates'])

        # add chromosome number, base pair distance, and p-value of the current chunk to results for all chunks
        self.append_to_results_all_chunks()

        # save results
        # save_process = multiprocessing.Process(target=self.save_results_regression)
        save_process = multiprocessing.Process(target=self.save_results_regression,
                                               args=(self.load('output_files')['output'][0],
                                                     self.load('config')['covariates'])
                                               )
        save_process.daemon = True
        save_process.start()
        save_process.join()
        save_process.terminate()

        # empty the dictionaries to release the memory because they are not needed anymore
        self.init_algorithm_attributes()

        # # if this is not the last chunk, set up the next chunk of SNPs
        # if not self.is_last_chunk():
        #     self.setup_next_chunk()
        # else:
        #     # if this is the last chunk, generate the manhattan plot first, and then, tell clients to download the results
        #     self.manhattan_plot()
        #     self.set_step(HyFedProjectStep.RESULT)

    def calculate_std_error_logistic_sub_chunk(self, start_index, end_index, queue_std_error):
        """ Compute logistic regression std error values for a sub-chunk """

        std_error_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.considered_snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_std_error.put(std_error_values)
                std_error_values = dict()

            if np.linalg.det(self.hessian_matrices[snp_index]) == 0:
                std_error_values[snp_index] = np.array(["NA" for _ in range(len(self.covariates) + 2)])
                continue

            std_error_values[snp_index] = np.sqrt(np.linalg.inv(self.hessian_matrices[snp_index]).diagonal())

        queue_std_error.put(std_error_values)
