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


@app_state('STD_Error_Linear', Role.BOTH)
class STDErrorLinear(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_STD_Error_Linear', Role.COORDINATOR)
        self.register_transition('Limbo', Role.PARTICIPANT)

    def run(self) -> str or None:
        load_attrs(self)
        self.beta_values = self.await_data()
        # self.sse_values = self.std_error_linear_step()
        self.std_error_linear_step()
        self.send_data_to_coordinator(data=self.sse_values, use_smpc=self.load('smpc_used'))
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_STD_Error_Linear'
        # It should not go to terminal until the last chunk
        return 'Limbo'

    # ##### functions related to the std error step of the linear regression algorithm
    def std_error_linear_step(self):
        """ Compute local sum square error values for the chunk """

        # set sse_values dictionary to empty at the beginning of the chunk
        self.sse_values = dict()

        # queue
        queue_sse = multiprocessing.Queue()

        # thread to read from the queues
        sse_read_thread = threading.Thread(target=self.read_queue_sse, args=(queue_sse,))
        sse_read_thread.daemon = True
        sse_read_thread.start()

        # global beta values
        # beta_values = self.global_parameters[SplinkGlobalParameter.BETA]
        self.current_chunk_size = len(self.beta_values)

        # processes to compute the local SSE values for the sub-chunks
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices,
                                                              self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_sse_values,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    self.beta_values, queue_sse))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        sse_read_thread.join()

        # close queue
        queue_sse.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionary to list
        sse_values = list()
        for snp_index in sorted(self.beta_values.keys()):
            sse_values.append(self.sse_values[snp_index])

        # python list of scalar values must be converted to a numpy array if compensator flag is set
        self.sse_values = np.array(sse_values)
        # return sse_values
        # share noisy local sse values with the server and noise with compensator
        # self.local_parameters[SplinkLocalParameter.SSE] = sse_values
        # self.set_compensator_flag({SplinkLocalParameter.SSE: DataType.NUMPY_ARRAY_FLOAT})

    def compute_sse_values(self, start_index, end_index, beta_values, queue_sse):
        """ Compute local sum square error value for a sub-chunk """

        sse_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in beta_values.keys():
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_sse.put(sse_values)
                sse_values = dict()

            # compute sum square error value for the SNP
            x_matrix, y_vector = self.get_x_matrix_y_vector(snp_index,
                                                            covariate=self.load('config')['covariates'],
                                                            algorithm=self.load('config')['algorithm'])
            beta_vector = beta_values[snp_index].reshape(-1, 1)
            y_predicted = np.dot(x_matrix, beta_vector)
            sse_values[snp_index] = np.sum(np.square(y_vector - y_predicted))

        queue_sse.put(sse_values)

    def read_queue_sse(self, queue_sse):
        while len(self.sse_values) < self.current_chunk_size:
            sse = queue_sse.get()
            self.sse_values.update(sse)


@app_state('Aggregate_STD_Error_Linear', Role.COORDINATOR)
class AggregateSTDErrorLinear(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Non_Missing_counts', Role.COORDINATOR)
        self.register_transition('terminal', Role.COORDINATOR)

    def run(self) -> str or None:
        load_attrs(self)
        sse_values = self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
        self.std_error_linear_step(sse_values)
        self.log("STD error linear step executed")
        if not self.is_last_chunk():
            data_to_send = self.setup_next_chunk()
            self.broadcast_data(data_to_send)
            share_attrs(self)
            return 'Non_Missing_counts'
        # if this is the last chunk, generate the manhattan plot first, and then, tell clients to download the results
        self.manhattan_plot(self.load('output_files')['manhattan'][0], self.log)
        # self.set_step(HyFedProjectStep.RESULT)
        return 'terminal'

    # ##### linear regression std-error step related functions
    def std_error_linear_step(self, sse_values):
        """ Compute linear regression standard error values using the aggregated SSE values """

        # aggregate SSE values from the clients
        # sse_values = self.compute_aggregated_parameter(SplinkLocalParameter.SSE, DataType.NUMPY_ARRAY_FLOAT)

        # convert sse list to dictionary
        self.sse_values = dict()
        snp_counter = -1
        for snp_index in sorted(self.considered_snp_indices):
            snp_counter += 1
            self.sse_values[snp_index] = sse_values[snp_counter]

        # initialize std_error_values as an empty dictionary
        self.std_error_values = dict()

        # queue
        queue_std_error = multiprocessing.Queue()

        # thread to read from the queue
        std_error_read_thread = threading.Thread(target=self.read_queue_std_error, args=(queue_std_error,))
        std_error_read_thread.daemon = True
        std_error_read_thread.start()
        self.log("read_queue_std_error executed")
        # processes to compute the std error values for the sub-chunks
        sub_chunk_start_indices, sub_chunk_end_indices = self.get_start_end_indices(cpu_cores=8)
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(sub_chunk_start_indices, sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.calculate_std_error_linear_sub_chunk,
                                              args=(start_index_sub_chunk,
                                                    end_index_sub_chunk,
                                                    queue_std_error,
                                                    self.load('config')['covariate_names']))
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

        # compute results (i.e. t-stats and p-values) for the chunk
        self.compute_results_regression(self.load('config')['algorithm'], self.load('config')['covariates'])

        # add chromosome number, base pair distance, and p-value of the current chunk to results for all chunks
        self.append_to_results_all_chunks()
        self.log(self.p_values[0])
        self.log(self.p_value_all_chunks[0])

        # save results
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

        # if this is not the last chunk, set up the next chunk of SNPs

    def calculate_std_error_linear_sub_chunk(self, start_index, end_index, queue_std_error, covariates):
        """ Compute linear regression std error values for a sub-chunk """

        std_error_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.considered_snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_std_error.put(std_error_values)
                std_error_values = dict()

            sigma_squared = self.sse_values[snp_index] / (self.non_missing_sample_counts[snp_index] - len(covariates) - 2)
            variance_values = (sigma_squared * self.xt_x_inverse_matrices[snp_index]).diagonal()
            std_error_values[snp_index] = np.sqrt(variance_values)

        queue_std_error.put(std_error_values)
