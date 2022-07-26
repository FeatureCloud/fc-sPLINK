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


@app_state('Beta_Linear', Role.BOTH)
class BetaLinear(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_Beta_Linear', Role.COORDINATOR)
        self.register_transition('STD_Error_Linear', Role.PARTICIPANT)

    def run(self) -> str or None:
        load_attrs(self)
        global_minor_allele_names, global_major_allele_names, self.snp_indices = self.await_data()
        self.minor_allele(global_minor_allele_names, global_major_allele_names)
        xt_x_matrix, xt_y_vector = self.beta_linear_step()
        self.send_data_to_coordinator(data=[xt_x_matrix, xt_y_vector], use_smpc=self.load('smpc_used'))
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_Beta_Linear'
        return 'STD_Error_Linear'

    def beta_linear_step(self):
        """ Compute X'X and X'Y matrices for the chunk """

        # set dictionaries to empty at the beginning of the chunk
        self.xt_x_matrix = dict()
        self.xt_y_vector = dict()

        # queues
        queue_xt_x_matrix = multiprocessing.Queue()
        queue_xt_y_vector = multiprocessing.Queue()

        # get SNP indices for which X'X and X'Y should be computed
        # self.snp_indices = self.global_parameters[SplinkGlobalParameter.SNP_INDEX]
        self.current_chunk_size = len(self.snp_indices)

        # threads to read from the queues
        xt_x_matrix_read_thread = threading.Thread(target=self.read_queue_xt_x_matrix,
                                                   args=(queue_xt_x_matrix,))
        xt_x_matrix_read_thread.daemon = True
        xt_x_matrix_read_thread.start()

        self.log(self.xt_x_matrix)

        xt_y_vector_read_thread = threading.Thread(target=self.read_queue_xt_y_vector,
                                                   args=(queue_xt_y_vector,))
        xt_y_vector_read_thread.daemon = True
        xt_y_vector_read_thread.start()

        # processes to compute local X'X and X'Y for the sub-chunks
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices, self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_beta_linear_parameters,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    queue_xt_x_matrix, queue_xt_y_vector))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read threads to be done
        xt_x_matrix_read_thread.join()
        xt_y_vector_read_thread.join()

        # close queues
        queue_xt_x_matrix.close()
        queue_xt_y_vector.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionaries to lists;
        xt_x_matrix = list()
        xt_y_vector = list()
        self.log(self.snp_indices)
        self.log(self.xt_x_matrix.keys())
        for snp_index in sorted(self.snp_indices):
            xt_x_matrix.append(self.xt_x_matrix[snp_index])
            xt_y_vector.append(self.xt_y_vector[snp_index])
        return xt_x_matrix, xt_y_vector
        # # share noisy local X'X matrix and X'Y vector with the server and noise with compensator
        # self.local_parameters[SplinkLocalParameter.XT_X_MATRIX] = xt_x_matrix
        # self.local_parameters[SplinkLocalParameter.XT_Y_VECTOR] = xt_y_vector
        # self.set_compensator_flag({SplinkLocalParameter.XT_X_MATRIX: DataType.LIST_NUMPY_ARRAY_FLOAT,
        #                            SplinkLocalParameter.XT_Y_VECTOR: DataType.LIST_NUMPY_ARRAY_FLOAT})

    def read_queue_xt_x_matrix(self, queue_xt_x_matrix):
        while len(self.xt_x_matrix) < self.current_chunk_size:
            xt_x = queue_xt_x_matrix.get()
            self.xt_x_matrix.update(xt_x)

    def read_queue_xt_y_vector(self, queue_xt_y_vector):
        while len(self.xt_y_vector) < self.current_chunk_size:
            xt_y = queue_xt_y_vector.get()
            self.xt_y_vector.update(xt_y)

    def compute_beta_linear_parameters(self, start_index, end_index, queue_xt_x_matrix, queue_xt_y_vector):
        """ Compute local X'X and X'Y for a sub-chunk """

        xt_x_matrices = dict()
        xt_y_vectors = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_xt_x_matrix.put(xt_x_matrices)
                queue_xt_y_vector.put(xt_y_vectors)
                xt_x_matrices = dict()
                xt_y_vectors = dict()
            x_matrix, y_vector = self.get_x_matrix_y_vector(snp_index,
                                                            covariate=self.load('config')['covariates'],
                                                            algorithm=self.load('config')['algorithm'])

            xt_x_matrices[snp_index] = np.dot(x_matrix.T, x_matrix)
            xt_y_vectors[snp_index] = np.dot(x_matrix.T, y_vector)

        queue_xt_x_matrix.put(xt_x_matrices)
        queue_xt_y_vector.put(xt_y_vectors)

    def minor_allele(self, global_minor_allele_names, global_major_allele_names):
        super().minor_allele(global_minor_allele_names, global_major_allele_names)
        self.attrs_to_share += ['snp_values', 'first_allele_names', 'second_allele_names']


@app_state('Aggregate_Beta_Linear', Role.COORDINATOR)
class AggregateBetaLinear(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('STD_Error_Linear', Role.COORDINATOR)

    def run(self) -> str or None:
        load_attrs(self)
        # xt_x_matrices, xt_y_vectors = self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
        if self.load('smpc_used'):
            xt_x_matrices, xt_y_vectors = self.await_data(n=1, unwrap=True, is_json=True)
        else:
            data = self.gather_data(is_json=False)
            xt_x_matrices, xt_y_vectors = _aggregate(data, SMPCOperation.ADD, n_parts=2)
        xt_x_matrices = np.array(xt_x_matrices) / len(self.clients)
        xt_y_vectors = np.array(xt_y_vectors) / len(self.clients)

        # beta_values = self.beta_linear_step(xt_x_matrices, xt_y_vectors)
        # self.broadcast_data(beta_values)
        self.beta_linear_step(xt_x_matrices, xt_y_vectors)
        self.broadcast_data(self.beta_values)
        self.attrs_to_share += ['considered_snp_indices', 'xt_x_matrices', 'xt_y_vectors', 'std_error_values','t_stat_values', 'p_values']
        share_attrs(self)
        return 'STD_Error_Linear'

    def beta_linear_step(self, xt_x_matrices, xt_y_vectors):
        """ Compute linear regression global beta values using the aggregated XTX and XTY matrices for the chunk """

        # aggregate X'X matrices and X'Y vectors from the clients
        # xt_x_matrices = self.compute_aggregated_parameter(SplinkLocalParameter.XT_X_MATRIX,
        #                                                   DataType.LIST_NUMPY_ARRAY_FLOAT)
        # xt_y_vectors = self.compute_aggregated_parameter(SplinkLocalParameter.XT_Y_VECTOR,
        #                                                  DataType.LIST_NUMPY_ARRAY_FLOAT)
        self.log("Beta Linear Step has Started")
        # convert lists to dictionaries
        self.xt_x_matrices = dict()
        self.xt_y_vectors = dict()
        snp_counter = -1
        for snp_index in sorted(self.considered_snp_indices.copy()):
            snp_counter += 1
            self.xt_x_matrices[snp_index] = xt_x_matrices[snp_counter]
            self.xt_y_vectors[snp_index] = xt_y_vectors[snp_counter]

        # initialize beta values and xt_x_inverse_matrices as empty dictionaries
        self.beta_values = dict()
        self.xt_x_inverse_matrices = dict()

        # queues
        queue_beta = multiprocessing.Queue()
        queue_xt_x_inverse = multiprocessing.Queue()

        # threads to read from the queue
        beta_read_thread = threading.Thread(target=self.read_queue_beta_linear, args=(queue_beta,))
        beta_read_thread.daemon = True
        beta_read_thread.start()

        xt_x_inverse_read_thread = threading.Thread(target=self.read_queue_xt_x_inverse, args=(queue_xt_x_inverse,))
        xt_x_inverse_read_thread.daemon = True
        xt_x_inverse_read_thread.start()

        # processes to compute the beta values and xt_x_inverse matrices for the sub-chunks
        sub_chunk_start_indices, sub_chunk_end_indices = self.get_start_end_indices(cpu_cores=8)
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(sub_chunk_start_indices, sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.calculate_beta_linear_sub_chunk,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    queue_beta, queue_xt_x_inverse,
                                                    self.load('config')['covariate_names']))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read threads to be done
        beta_read_thread.join()
        xt_x_inverse_read_thread.join()

        # close queues
        queue_beta.close()
        queue_xt_x_inverse.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # update considered index set
        for snp_index in self.considered_snp_indices:
            if self.beta_values[snp_index][0] == "NA":
                self.considered_snp_indices.discard(snp_index)
                self.std_error_values[snp_index] = self.beta_values[snp_index]
                self.t_stat_values[snp_index] = self.beta_values[snp_index]
                self.p_values[snp_index] = self.beta_values[snp_index]
                continue

        # only share beta values for considered SNPs with clients to compute sum square error values
        self.beta_values = {snp_index: self.beta_values[snp_index] for snp_index in self.considered_snp_indices}
        # return beta_values
        # self.global_parameters[SplinkGlobalParameter.BETA] = beta_values
        #
        # # tell clients to go to std-error step
        # self.set_step(SplinkProjectStep.STD_ERROR_LINEAR)

    def read_queue_xt_x_inverse(self, queue_xt_x_inverse):
        while len(self.xt_x_inverse_matrices) < len(self.considered_snp_indices):
            xt_x_inverse = queue_xt_x_inverse.get()
            self.xt_x_inverse_matrices.update(xt_x_inverse)

    def read_queue_beta_linear(self, queue_beta_linear):
        while len(self.beta_values) < len(self.considered_snp_indices):
            betas = queue_beta_linear.get()
            self.beta_values.update(betas)

    def calculate_beta_linear_sub_chunk(self, start_index, end_index, queue_beta, queue_xt_x_inverse, covariates):
        """ Compute linear regression beta values for a sub-chunk """

        beta_values = dict()
        xt_x_inverse_matrices = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.considered_snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_beta.put(beta_values)
                queue_xt_x_inverse.put(xt_x_inverse_matrices)
                beta_values = dict()
                xt_x_inverse_matrices = dict()

            if np.linalg.det(self.xt_x_matrices[snp_index]) == 0:
                self.beta_values[snp_index] = np.array(["NA" for _ in range(len(covariates) + 2)])
                continue

            xt_x_inverse_matrices[snp_index] = np.linalg.inv(self.xt_x_matrices[snp_index])
            beta_values[snp_index] = np.dot(xt_x_inverse_matrices[snp_index], self.xt_y_vectors[snp_index]).flatten()

        queue_beta.put(beta_values)
        queue_xt_x_inverse.put(xt_x_inverse_matrices)


def _aggregate(data, operation: SMPCOperation, n_parts: int = 1):
    """
    Aggregates a list of received values.

    Parameters
    ----------
    data : array_like
        list of data pieces
    operation : SMPCOperation
        operation to use for aggregation (add or multiply)

    Returns
    ----------
    aggregated value
    """

    def aggregate(data):
        data_np = [np.array(d) for d in data]

        aggregate = data_np[0]

        if operation == SMPCOperation.ADD:
            for d in data_np[1:]:
                aggregate = aggregate + d

        if operation == SMPCOperation.MULTIPLY:
            for d in data_np[1:]:
                aggregate = aggregate * d
        return aggregate

    if n_parts == 1:
        return aggregate(data)
    return [aggregate([data[i][n] for i in range(len(data))]) for n in range(n_parts)]
