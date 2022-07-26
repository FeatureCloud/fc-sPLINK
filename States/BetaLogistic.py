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


@app_state('Beta_Logistic', Role.BOTH)
class BetaLogistic(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_Beta_Logistic', Role.COORDINATOR)
        # Self loop
        self.register_transition('Beta_Logistic', Role.PARTICIPANT)
        self.register_transition('STD_Error_Logistic', Role.BOTH)

    def run(self) -> str or None:
        load_attrs(self)
        # beta_values, self.current_beta_iteration = self.await_data()
        received_data = self.await_data()
        if len(received_data) == 2:
            beta_values, converged = received_data
            self.beta_values = beta_values
            if converged:
                self.store('beta_values', beta_values)
                return 'STD_Error_Logistic'
        else:
            global_minor_allele_names, global_major_allele_names, beta_values, converged = received_data
            self.beta_values = beta_values
            self.minor_allele(global_minor_allele_names, global_major_allele_names)


        data_to_send = self.beta_logistic_step(beta_values)
        self.send_data_to_coordinator(data_to_send)
        share_attrs(self)
        if self.is_coordinator:
            return 'Aggregate_Beta_Logistic'
        return 'Beta_Logistic'

    # ##### logistic regression beta step related functions
    def beta_logistic_step(self, beta_values):
        """ Compute local gradient and Hessian matrices as well as log likelihood values for the chunk """

        # set gradient, Hessian, and log likelihood dictionaries to empty at the beginning of the chunk
        self.gradient_vectors = dict()
        self.hessian_matrices = dict()
        self.log_likelihood_values = dict()

        # queues
        queue_gradient = multiprocessing.Queue()
        queue_hessian = multiprocessing.Queue()
        queue_log_likelihood = multiprocessing.Queue()

        # thread to read from the queues
        gradient_read_thread = threading.Thread(target=self.read_queue_gradient, args=(queue_gradient,))
        gradient_read_thread.daemon = True
        gradient_read_thread.start()

        hessian_read_thread = threading.Thread(target=self.read_queue_hessian, args=(queue_hessian,))
        hessian_read_thread.daemon = True
        hessian_read_thread.start()

        log_likelihood_read_thread = threading.Thread(target=self.read_queue_log_likelihood,
                                                      args=(queue_log_likelihood,))
        log_likelihood_read_thread.daemon = True
        log_likelihood_read_thread.start()

        # global beta values and current beta iteration
        self.current_chunk_size = len(beta_values)

        # processes to compute the gradient, Hessian, and log likelihood values for sub-chunks
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(self.sub_chunk_start_indices,
                                                              self.sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.compute_beta_logistic_parameters,
                                              args=(start_index_sub_chunk, end_index_sub_chunk,
                                                    beta_values, queue_gradient, queue_hessian,
                                                    queue_log_likelihood))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        gradient_read_thread.join()
        hessian_read_thread.join()
        log_likelihood_read_thread.join()

        # close queues
        queue_gradient.close()
        queue_hessian.close()
        queue_log_likelihood.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # convert dictionary to list
        gradient_vectors = list()
        hessian_matrices = list()
        log_likelihood_values = list()
        for snp_index in sorted(beta_values.keys()):
            gradient_vectors.append(self.gradient_vectors[snp_index])
            hessian_matrices.append(self.hessian_matrices[snp_index])
            log_likelihood_values.append(self.log_likelihood_values[snp_index])

        # python list of scalars must be converted to a numpy array if compensator flag is set
        log_likelihood_values = np.array(log_likelihood_values)
        return gradient_vectors, hessian_matrices, log_likelihood_values
        # # share the noisy local gradient, Hessian, and log likelihood values with the server and noise with compensator
        # self.local_parameters[SplinkLocalParameter.GRADIENT] = gradient_vectors
        # self.local_parameters[SplinkLocalParameter.HESSIAN] = hessian_matrices
        # self.local_parameters[SplinkLocalParameter.LOG_LIKELIHOOD] = log_likelihood_values
        # self.set_compensator_flag({SplinkLocalParameter.GRADIENT: DataType.LIST_NUMPY_ARRAY_FLOAT,
        #                            SplinkLocalParameter.HESSIAN: DataType.LIST_NUMPY_ARRAY_FLOAT,
        #                            SplinkLocalParameter.LOG_LIKELIHOOD: DataType.NUMPY_ARRAY_FLOAT})

    def compute_beta_logistic_parameters(self, start_index, end_index, beta_values, queue_gradient, queue_hessian,
                                         queue_log_likelihood):
        """ Compute local gradient vector, Hessian matrix, and log likelihood values for a sub-chunk """

        epsilon = np.finfo(float).eps  # to avoid log(0)
        gradient_vectors = dict()
        hessian_matrices = dict()
        log_likelihood_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in beta_values.keys():
                continue

            # put results in the queues whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_gradient.put(gradient_vectors)
                queue_hessian.put(hessian_matrices)
                queue_log_likelihood.put(log_likelihood_values)
                gradient_vectors = dict()
                hessian_matrices = dict()
                log_likelihood_values = dict()

            x_matrix, y_vector = self.get_x_matrix_y_vector(snp_index,
                                                            covariate=self.load('config')['covariates'],
                                                            algorithm=self.load('config')['algorithm'])
            beta_vector = beta_values[snp_index].reshape(-1, 1)

            # gradient
            x_beta_product = np.dot(x_matrix, beta_vector)
            y_predicted = 1 / (1 + np.exp(-x_beta_product))
            gradient_vectors[snp_index] = np.dot(x_matrix.T, (y_vector - y_predicted))

            # hessian matrix
            hessian_matrices[snp_index] = np.dot(np.multiply(x_matrix.T, (y_predicted * (1 - y_predicted)).T),
                                                 x_matrix)

            # log likelihood
            log_likelihood_values[snp_index] = np.sum(
                y_vector * np.log(y_predicted + epsilon) + (1 - y_vector) * np.log(1 - y_predicted + epsilon))

        queue_gradient.put(gradient_vectors)
        queue_hessian.put(hessian_matrices)
        queue_log_likelihood.put(log_likelihood_values)

    def minor_allele(self, global_minor_allele_names, global_major_allele_names):
        super().minor_allele(global_minor_allele_names, global_major_allele_names)
        self.attrs_to_share += ['snp_values', 'first_allele_names', 'second_allele_names']


@app_state('Aggregate_Beta_Logistic', Role.COORDINATOR)
class BetaLogistic(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Beta_Logistic', Role.COORDINATOR)

    def run(self) -> str or None:
        # gradient_vectors, hessian_matrices, log_likelihood_values = self.gather_data()
        gradient_vectors, hessian_matrices, log_likelihood_values = \
            self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))
        load_attrs(self)
        self.beta_logistic_step(gradient_vectors, hessian_matrices, log_likelihood_values)
        is_converged = self.is_converged()
        if is_converged:
            self.current_beta_iteration += 1
            beta_values = {snp_index: self.beta_values[snp_index] for snp_index in
                           self.considered_in_process_snp_indices}

        else:
            beta_values = {snp_index: self.beta_values[snp_index] for snp_index in self.considered_snp_indices}
        self.broadcast_data(data=[beta_values, is_converged])
        share_attrs(self)
        return 'Beta_Logistic'

    def is_converged(self):
        return self.current_beta_iteration != self.load('config')['max_iterations'] and len(
            self.considered_in_process_snp_indices) != 0

    def beta_logistic_step(self, gradient_vectors, hessian_matrices, log_likelihood_values):
        """ Compute logistic regression global beta values using the aggregated gradient and Hessian matrices for the chunk """

        # convert lists to dictionaries
        self.gradient_vectors = dict()
        self.hessian_matrices = dict()
        self.new_log_likelihood_values = dict()
        snp_counter = -1
        for snp_index in sorted(self.considered_in_process_snp_indices):
            snp_counter += 1
            self.gradient_vectors[snp_index] = gradient_vectors[snp_counter]
            self.hessian_matrices[snp_index] = hessian_matrices[snp_counter]
            self.new_log_likelihood_values[snp_index] = log_likelihood_values[snp_counter]

        # initialize new beta values as an empty dictionary
        self.new_beta_values = dict()

        # queue
        queue_beta_values = multiprocessing.Queue()

        # thread to read from the queue
        beta_value_read_thread = threading.Thread(target=self.read_queue_beta_logistic, args=(queue_beta_values,))
        beta_value_read_thread.daemon = True
        beta_value_read_thread.start()

        # processes to compute the new beta values for the sub-chunks
        sub_chunk_start_indices, sub_chunk_end_indices = self.get_start_end_indices(cpu_cores=8)
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(sub_chunk_start_indices, sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.calculate_beta_logistic_sub_chunk,
                                              args=(start_index_sub_chunk, end_index_sub_chunk, queue_beta_values,))
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        beta_value_read_thread.join()

        # close queues
        queue_beta_values.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()

        # update beta values
        for snp_index in self.new_beta_values.keys():
            self.beta_values[snp_index] = self.new_beta_values[snp_index]

        # update considered index set
        for snp_index in self.considered_in_process_snp_indices:
            if self.beta_values[snp_index][0] == "NA":
                self.considered_snp_indices.discard(snp_index)
                self.std_error_values[snp_index] = self.beta_values[snp_index]
                self.t_stat_values[snp_index] = self.beta_values[snp_index]
                self.p_values[snp_index] = self.beta_values[snp_index]
                continue

        # check whether beta values for the SNP converged. If so, remove the SNP index from the in_process indices
        for snp_index in self.considered_in_process_snp_indices:
            old_log_likelihood = self.log_likelihood_values[snp_index]
            new_log_likelihood = self.new_log_likelihood_values[snp_index]
            if self.has_converged(old_log_likelihood, new_log_likelihood):
                self.in_process_snp_indices.discard(snp_index)

        # update log likelihood values
        for snp_index in self.new_log_likelihood_values.keys():
            self.log_likelihood_values[snp_index] = self.new_log_likelihood_values[snp_index]

        # if there are still SNPs whose beta values not converged and max iterations not reached yet,
        # share updated global beta values (excluding those ignored or converged) with the clients and stay in beta_logistic step
        self.considered_in_process_snp_indices = self.considered_snp_indices.intersection(self.in_process_snp_indices)
        # if self.current_beta_iteration != self.max_iterations and len(self.considered_in_process_snp_indices) != 0:
        #     self.current_beta_iteration += 1
        #     beta_values = {snp_index: self.beta_values[snp_index] for snp_index in
        #                    self.considered_in_process_snp_indices}
        #     return beta_values, self.current_beta_iteration
        #     # self.global_parameters[SplinkGlobalParameter.BETA] = beta_values
        #     # self.global_parameters[SplinkGlobalParameter.CURRENT_BETA_ITERATION] = self.current_beta_iteration
        #     # logger.debug(f'Project {self.project_id}: Beta iteration # {self.current_beta_iteration} done!')
        #
        # # if beta max iterations reached or all beta values converged, share updated beta values (excluding ignored SNPs)
        # # with clients and go to the std_error_logistic step
        # beta_values = {snp_index: self.beta_values[snp_index] for snp_index in self.considered_snp_indices}
        # return beta_values
        # # else:
        # #     beta_values = {snp_index: self.beta_values[snp_index] for snp_index in self.considered_snp_indices}
        # #     self.global_parameters[SplinkGlobalParameter.BETA] = beta_values
        # #     self.set_step(SplinkProjectStep.STD_ERROR_LOGISTIC)

    def calculate_beta_logistic_sub_chunk(self, start_index, end_index, queue_beta_values):
        """ Compute logistic regression beta values for a sub-chunk """

        new_beta_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.considered_in_process_snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_beta_values.put(new_beta_values)
                new_beta_values = dict()

            if np.linalg.det(self.hessian_matrices[snp_index]) == 0:
                new_beta_values[snp_index] = np.array(["NA" for _ in range(len(self.covariates) + 2)])
                continue

            hessian_inverse_matrix = np.linalg.inv(self.hessian_matrices[snp_index])
            beta_update_vector = np.dot(hessian_inverse_matrix, self.gradient_vectors[snp_index])
            new_beta_vector = self.beta_values[snp_index].reshape(-1, 1) + beta_update_vector
            new_beta_values[snp_index] = new_beta_vector.flatten()

        queue_beta_values.put(new_beta_values)

    def read_queue_beta_logistic(self, queue_beta_values):
        while len(self.new_beta_values) < len(self.considered_in_process_snp_indices):
            new_betas = queue_beta_values.get()
            self.new_beta_values.update(new_betas)

    def has_converged(self, old_log_likelihood, new_log_likelihood):
        """ Determine whether beta values has converged based on the old and new values of log likelihood """

        if old_log_likelihood is None:
            return False

        delta_log_likelihood = np.abs(old_log_likelihood - new_log_likelihood)
        if delta_log_likelihood > self.delta_log_likelihood_threshold:
            return False

        return True
