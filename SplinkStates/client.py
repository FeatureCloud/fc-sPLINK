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
