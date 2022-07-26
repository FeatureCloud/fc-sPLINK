from SplinkStates.splink import Splink
import numpy as np
from States import ALGORITHM
from Utils.gwas_dataset import MissingValue


class SplinkClient(Splink):
    def __init__(self):
        super().__init__()

    def get_x_matrix_y_vector(self, snp_index, covariate, algorithm):
        """ Create feature matrix and label vector """
        # get non-missing rows after considering SNP values
        snp_indices_non_missing = self.snp_values[snp_index] != MissingValue.SNP
        index_values_non_missing = np.logical_and(self.non_missing_index_values, snp_indices_non_missing)

        # create feature matrix
        snp_vector = self.snp_values[snp_index][index_values_non_missing].reshape(-1, 1).astype(np.int64)

        if len(covariate) == 0:
            x_matrix = np.concatenate((np.ones((len(snp_vector), 1)).astype(np.int64), snp_vector), axis=1)
        else:
            x_matrix = np.concatenate((np.ones((len(snp_vector), 1)).astype(np.int64),
                                       snp_vector, self.covariate_matrix[index_values_non_missing, :]), axis=1)

        # create label vector
        y_vector = self.phenotype_values[index_values_non_missing].reshape(-1, 1)

        if algorithm in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]:
            y_vector = y_vector.astype(np.uint8)

        return x_matrix, y_vector

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

    # ##### Queue functions
    def read_queue_non_missing(self, queue_non_missing):
        while len(self.non_missing_sample_counts) < self.current_chunk_size:
            sample_count_non_missing = queue_non_missing.get()
            self.non_missing_sample_counts.update(sample_count_non_missing)

    def read_queue_allele_counts(self, queue_allele_counts):
        while len(self.allele_counts) < self.current_chunk_size:
            count_alleles = queue_allele_counts.get()
            self.allele_counts.update(count_alleles)

    def read_queue_contingency_tables(self, queue_contingency_tables):
        while len(self.contingency_tables) < self.current_chunk_size:
            cont_table = queue_contingency_tables.get()
            self.contingency_tables.update(cont_table)

    def read_queue_xt_x_matrix(self, queue_xt_x_matrix):
        while len(self.xt_x_matrix) < self.current_chunk_size:
            xt_x = queue_xt_x_matrix.get()
            self.xt_x_matrix.update(xt_x)

    def read_queue_xt_y_vector(self, queue_xt_y_vector):
        while len(self.xt_y_vector) < self.current_chunk_size:
            xt_y = queue_xt_y_vector.get()
            self.xt_y_vector.update(xt_y)

    def read_queue_sse(self, queue_sse):
        while len(self.sse_values) < self.current_chunk_size:
            sse = queue_sse.get()
            self.sse_values.update(sse)

    def read_queue_gradient(self, queue_gradient):
        while len(self.gradient_vectors) < self.current_chunk_size:
            gradient = queue_gradient.get()
            self.gradient_vectors.update(gradient)

    def read_queue_hessian(self, queue_hessian):
        while len(self.hessian_matrices) < self.current_chunk_size:
            hessian_matrix = queue_hessian.get()
            self.hessian_matrices.update(hessian_matrix)

    def read_queue_log_likelihood(self, queue_log_likelihood):
        while len(self.log_likelihood_values) < self.current_chunk_size:
            log_likelihood = queue_log_likelihood.get()
            self.log_likelihood_values.update(log_likelihood)