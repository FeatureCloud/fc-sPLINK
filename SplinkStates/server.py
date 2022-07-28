import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Utils.gwas_dataset import MissingValue
from States import ALGORITHM
import multiprocessing
import threading
from scipy import stats
from SplinkStates.splink import Splink


class SplinkServer(Splink):
    def __init__(self):
        super().__init__()

    # ##### Helper functions

    def compute_p_values(self, algorithm, covariate):
        """ Compute p-values for a chunk with multi-processing """
        queue_p_values = multiprocessing.Queue()

        # thread to read from the queue
        p_value_read_thread = threading.Thread(target=self.read_queue_p_values, args=(queue_p_values,))
        p_value_read_thread.daemon = True
        p_value_read_thread.start()

        # processes to compute the p-values for the sub-chunks
        sub_chunk_start_indices, sub_chunk_end_indices = self.get_start_end_indices(cpu_cores=8)
        process_list = list()
        for start_index_sub_chunk, end_index_sub_chunk in zip(sub_chunk_start_indices, sub_chunk_end_indices):
            process = multiprocessing.Process(target=self.calculate_p_values_sub_chunk,
                                              args=(start_index_sub_chunk, end_index_sub_chunk, queue_p_values,
                                                    algorithm, covariate)
                                              )
            process_list.append(process)
            process.daemon = True
            process.start()

        # wait for read thread to be done
        p_value_read_thread.join()

        # close queues
        queue_p_values.close()

        # terminate the processes
        for proc in process_list:
            proc.terminate()
        return f"p-value computation is done for chunk # {self.current_chunk}!"

    def calculate_p_values_sub_chunk(self, start_index, end_index, queue_p_values, algorithm, covariate):
        """ Compute p-values for a sub-chunk """

        p_values = dict()

        for snp_index in np.arange(start_index, end_index):
            if snp_index not in self.considered_snp_indices:
                continue

            # put results in the queue whenever computation is done for 1000 SNPs
            if snp_index % 1001 == 1000:
                queue_p_values.put(p_values)
                p_values = dict()

            if algorithm == ALGORITHM.CHI_SQUARE:
                p_values[snp_index] = 1 - stats.chi2.cdf(self.chi_square_values[snp_index], 1)
            elif algorithm == ALGORITHM.LINEAR_REGRESSION:
                degree_of_freedom = self.non_missing_sample_counts[snp_index] - len(
                    covariate) - 2
                p_values[snp_index] = 2 * (1 - stats.t.cdf(np.abs(self.t_stat_values[snp_index]), degree_of_freedom))
            elif algorithm == ALGORITHM.LOGISTIC_REGRESSION:
                p_values[snp_index] = 1 - stats.chi2.cdf(np.square(np.array(self.t_stat_values[snp_index])), 1)

        queue_p_values.put(p_values)

    def read_queue_p_values(self, queue_p_values):
        while len(self.p_values) < len(self.considered_snp_indices):
            prob_values = queue_p_values.get()
            self.p_values.update(prob_values)

    #
    # # ##### Chi-square result computation/saving functions
    def compute_maf(self):
        """ Compute minor allele frequency of case/control for the chunk """
        for snp_index in self.considered_snp_indices:
            minor_case = self.contingency_tables[snp_index][0]
            major_case = self.contingency_tables[snp_index][1]
            minor_control = self.contingency_tables[snp_index][2]
            major_control = self.contingency_tables[snp_index][3]

            self.maf_case[snp_index] = minor_case / (minor_case + major_case)
            self.maf_control[snp_index] = minor_control / (minor_control + major_control)
        self.attrs_to_share += ['maf_case', 'maf_control']
        return f'case/control minor allele frequency computation is done for chunk # {self.current_chunk}!'

    def compute_chi_square_values(self):
        """ Compute chi-square value for the chunk """
        for snp_index in self.considered_snp_indices:
            # observed allele counts
            observed_allele_counts = self.contingency_tables[snp_index]

            # expected allele counts
            expected_allele_counts = np.zeros(4)
            case_count = self.contingency_tables[snp_index][0] + self.contingency_tables[snp_index][1]
            control_count = self.contingency_tables[snp_index][2] + self.contingency_tables[snp_index][3]
            minor_count = self.contingency_tables[snp_index][0] + self.contingency_tables[snp_index][2]
            major_count = self.contingency_tables[snp_index][1] + self.contingency_tables[snp_index][3]
            total_count = case_count + control_count

            expected_allele_counts[0] = (case_count * minor_count) / total_count
            expected_allele_counts[1] = (case_count * major_count) / total_count
            expected_allele_counts[2] = (control_count * minor_count) / total_count
            expected_allele_counts[3] = (control_count * major_count) / total_count

            # compute chi-square value
            chi_square = np.sum(np.square(observed_allele_counts - expected_allele_counts) / expected_allele_counts)
            self.chi_square_values[snp_index] = chi_square
        self.attrs_to_share += ['chi_square_values']

        return f"chi-square computation is done for chunk # {self.current_chunk}!"

    def compute_odd_ratio_values(self):
        """ Compute odd ratio value for the chunk """
        for snp_index in self.considered_snp_indices:
            minor_case = self.contingency_tables[snp_index][0]
            major_case = self.contingency_tables[snp_index][1]
            minor_control = self.contingency_tables[snp_index][2]
            major_control = self.contingency_tables[snp_index][3]

            if (major_case * minor_control) != 0:
                self.odd_ratio_values[snp_index] = (minor_case * major_control) / (major_case * minor_control)
            else:
                self.odd_ratio_values[snp_index] = "NA"
        self.attrs_to_share += ['odd_ratio_values']
        return f"odd-ratio computation is done for chunk # {self.current_chunk}!"

    #
    # def compute_results_chi_square(self):
    #     """ Compute MAF for case/control, chi-square, odd-ratio, and p-values for chi-square algorithm """
    #
    #     try:
    #         self.compute_maf()
    #         self.compute_chi_square_values()
    #         self.compute_odd_ratio_values()
    #         self.compute_p_values()
    #     except Exception as result_computation_error:
    #         logger.error(f"Chi-square result computation error: {result_computation_error}")
    #         self.project_failed()
    #
    def save_results_chi_square(self, output_file):
        """ Save chi-square algorithm results for the chunk into the file """

        # create result directory/file if they do not already exist
        result_file = open(output_file, 'a')

        # write the result file header in the first chunk
        if self.current_chunk == 1:
            result_file.write('CHR,SNP,BP,A1,F_A,F_U,A2,CHISQ,P,OR')

        for snp_index in np.arange(self.chunk_start_index, self.chunk_end_index):
            snp_id = self.snp_id_values[snp_index].decode('utf-8')
            chromosome_number, snp_name, base_pair_distance = snp_id.split('\t')
            minor_allele = self.minor_allele_names[snp_index]
            major_allele = self.major_allele_names[snp_index]
            maf_case = round_result(self.maf_case[snp_index])
            maf_control = round_result(self.maf_control[snp_index])
            chi_square = round_result(self.chi_square_values[snp_index])
            p_value = round_result(self.p_values[snp_index])
            odd_ratio = round_result(self.odd_ratio_values[snp_index])

            csv_row = f'{chromosome_number},{snp_name},{base_pair_distance},{minor_allele},{maf_case},' \
                      f'{maf_control},{major_allele},{chi_square},{p_value},{odd_ratio}'

            result_file.write("\n" + str(csv_row))

        result_file.close()

        print(f'Saving results done for chunk # {self.current_chunk}!')

    #
    # # ###### Linear/logistic regression result computation/saving functions
    def compute_t_stat_values(self):
        """ Compute T statistics for the chunk """

        for snp_index in self.considered_snp_indices:
            self.t_stat_values[snp_index] = self.beta_values[snp_index] / self.std_error_values[snp_index]

        print(f'T statistics computation done for chunk # {self.current_chunk}!')

    def compute_results_regression(self, algorithm, covariates):
        """ Compute t-stat and p-values for the linear/logistic regression algorithm """
        self.compute_t_stat_values()
        self.compute_p_values(algorithm, covariates)

    # # ############## Chunking functions
    # def init_chunks(self):
    #     """ Set the total number of chunks and start/end indices of the chunks """
    #
    #     try:
    #         self.total_chunks = int(np.ceil(len(self.snp_id_values) / self.chunk_size))
    #         for split in np.array_split(np.arange(len(self.snp_id_values)), self.total_chunks):
    #             self.start_indices_chunks.append(split[0])
    #             self.end_indices_chunks.append(split[-1] + 1)
    #
    #         logger.debug(f'Project {self.project_id}: Initializing of chunks is done!')
    #
    #     except Exception as init_chunk_exp:
    #         logger.error(f'Project {self.project_id}: {init_chunk_exp}')
    #         self.project_failed()
    #
    def init_algorithm_attributes(self):
        """ Set the chi-square or linear/logistic regression algorithm related dictionaries to empty """

        self.non_missing_sample_counts = dict()
        self.allele_counts = dict()
        self.minor_allele_names = dict()
        self.major_allele_names = dict()
        self.minor_allele_counts = dict()
        self.major_allele_counts = dict()
        self.minor_allele_frequencies = dict()
        self.major_allele_frequencies = dict()

        self.contingency_tables = dict()
        self.maf_case = dict()
        self.maf_control = dict()
        self.chi_square_values = dict()
        self.odd_ratio_values = dict()

        self.xt_x_matrices = dict()
        self.xt_y_vectors = dict()
        self.xt_x_inverse_matrices = dict()
        self.sse_values = dict()

        self.gradient_vectors = dict()
        self.hessian_matrices = dict()
        self.new_log_likelihood_values = dict()
        self.new_beta_values = dict()
        self.log_likelihood_values = dict()

        self.beta_values = dict()
        self.std_error_values = dict()
        self.t_stat_values = dict()
        self.p_values = dict()

    def setup_next_chunk(self):
        """ For the next chunk of SNPs:
                set the start/end chunk index, increment chunk number,
                set the chunk related global parameter values, and go to non-missing-count step
        """
        # set the chunk attribute values
        self.chunk_start_index = self.start_indices_chunks[self.current_chunk]
        self.chunk_end_index = self.end_indices_chunks[self.current_chunk]
        self.current_chunk += 1
        self.considered_snp_indices = set(np.arange(self.chunk_start_index, self.chunk_end_index)).copy()
        self.in_process_snp_indices = set(
            np.arange(self.chunk_start_index, self.chunk_end_index)).copy()  # used in BETA step of logistic reg
        self.current_beta_iteration = 1

        print(f'Chunk # {self.current_chunk} initialized!')
        data_to_send = [self.current_chunk, self.total_chunks, self.considered_snp_indices, self.chunk_start_index,
                        self.chunk_end_index]
        return data_to_send

    # ############## Helper functions
    def is_last_chunk(self):
        """ Check whether current chunk is the last one """

        return self.current_chunk == self.total_chunks

    #
    # def has_converged(self, old_log_likelihood, new_log_likelihood):
    #     """ Determine whether beta values has converged based on the old and new values of log likelihood """
    #
    #     try:
    #         if old_log_likelihood is None:
    #             return False
    #
    #         delta_log_likelihood = np.abs(old_log_likelihood - new_log_likelihood)
    #         if delta_log_likelihood > self.delta_log_likelihood_threshold:
    #             return False
    #
    #         return True
    #     except Exception as convergence_exception:
    #         logger.error(f'Project {self.project_id}: {convergence_exception}')
    #         self.project_failed()
    #
    def get_start_end_indices(self, cpu_cores):
        """ Determine start/end indices for sub-chunks assigned to each process/core """
        chunk_size = self.chunk_end_index - self.chunk_start_index

        # ensure each process/core will compute at least one SNP statistics
        if chunk_size < cpu_cores:
            cpu_cores = 1

        sub_chunk_size = int(np.ceil(chunk_size / cpu_cores))
        start_indices = np.arange(self.chunk_start_index, self.chunk_end_index, sub_chunk_size)
        end_indices = start_indices + sub_chunk_size
        end_indices[-1] = self.chunk_end_index

        return start_indices, end_indices

    def append_to_results_all_chunks(self):
        """ Add the chromosome numbers, base pair distances, and p-values of the current chunk to
            the corresponding lists for all chunks """

        for snp_index in np.arange(self.chunk_start_index, self.chunk_end_index):
            snp_id = self.snp_id_values[snp_index].decode('utf-8')
            chromosome_number, _, base_pair_distance = snp_id.split('\t')
            p_value = round_result(self.p_values[snp_index])

            self.chromosome_number_all_chunks.append(chromosome_number)
            self.base_pair_distance_all_chunks.append(base_pair_distance)
            self.p_value_all_chunks.append(p_value)
        self.attrs_to_share += ['chromosome_number_all_chunks', 'base_pair_distance_all_chunks', 'p_value_all_chunks']

    def manhattan_plot(self, manhatan_filename, log):
        """ draw Manhattan plot for p-values after processing of all chunks finished """
        manhattan_dict = {'CHR': self.chromosome_number_all_chunks,
                          'BP': self.base_pair_distance_all_chunks,
                          'P': self.p_value_all_chunks}
        log(manhattan_dict['CHR'])
        log(manhattan_dict['BP'])
        log(manhattan_dict['P'])
        manhattan_df = pd.DataFrame.from_dict(manhattan_dict)
        manhattan_df['P_LOG10'] = [-np.log10(row) for row in manhattan_df.P.values]

        manhattan_df.CHR = manhattan_df.CHR.astype('category')
        manhattan_df.CHR = manhattan_df.CHR.cat.set_categories(list(set(manhattan_df.CHR)), ordered=True)
        manhattan_df = manhattan_df.sort_values(['CHR', 'BP'])
        manhattan_df['ind'] = range(len(manhattan_df))
        manhattan_df_grouped = manhattan_df.groupby('CHR')
        fig = plt.figure(figsize=(24, 8), dpi=80)
        ax = fig.add_subplot(111)
        colors = ['blue', 'green', 'purple', 'brown']
        x_labels = []
        x_labels_pos = []
        try:
            for num, (name, group) in enumerate(manhattan_df_grouped):
                group.plot(kind='scatter', x='ind', y='P_LOG10', color=colors[num % len(colors)], ax=ax)

                x_labels.append(name)
                x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

                print(x_labels_pos[-1])
            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels)
            ax.set_xlim([0, len(manhattan_df)])
            ax.set_xlabel('Chromosome')
            ax.set_ylabel('-log10(p)')

            plt.savefig(manhatan_filename, format='png')

            log(f'Manhattan plot created!')
        except Exception as m_exp:
            log("Manhattan plot is not supported for this scenario")

    # used in std-error step of linear/logistic regression
    def read_queue_std_error(self, queue_std_error):
        while len(self.std_error_values) < len(self.considered_snp_indices):
            std_error = queue_std_error.get()
            self.std_error_values.update(std_error)

    def save_results_regression(self, res_filename, covariates):
        """ Save the linear/logistic regression results for the chunk into the file """

        # create result directory/file if they do not already exist
        # result_dir = self.create_result_dir()

        # if algorithm == ALGORITHM.LINEAR_REGRESSION:
        #     result_file = open(f'{result_dir}/linear-regression-result.csv', 'a')
        # else:
        #     result_file = open(f'{result_dir}/logistic-regression-result.csv', 'a')

        result_file = open(res_filename, 'a')

        # write the result file header in the first chunk
        if self.current_chunk == 1:
            result_file.write('CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P')
        for snp_index in np.arange(self.chunk_start_index, self.chunk_end_index):
            snp_id = self.snp_id_values[snp_index].decode('utf-8')
            chromosome_number, snp_name, base_pair_distance = snp_id.split('\t')

            beta_counter = 1
            minor_allele = self.minor_allele_names[snp_index]
            feature_name = 'ADD'
            non_missing_samples = round_result(self.non_missing_sample_counts[snp_index])

            beta_value = round_result(self.beta_values[snp_index][beta_counter])
            t_stat_value = round_result(self.t_stat_values[snp_index][beta_counter])
            p_value = round_result(self.p_values[snp_index][beta_counter])

            csv_row = f'{chromosome_number},{snp_name},{base_pair_distance},' \
                      f'{minor_allele},{feature_name},{non_missing_samples},' \
                      f'{beta_value},{t_stat_value},{p_value}'

            result_file.write("\n" + str(csv_row))

            for covariate in covariates:
                # for beta_counter, covariate in enumerate(covariates):
                beta_counter += 1
                beta_value = round_result(self.beta_values[snp_index][beta_counter])
                t_stat_value = round_result(self.t_stat_values[snp_index][beta_counter])
                p_value = round_result(self.p_values[snp_index][beta_counter])

                csv_row = f'{chromosome_number},{snp_name},{base_pair_distance},' \
                          f'{minor_allele},{covariate},{non_missing_samples},' \
                          f'{beta_value},{t_stat_value},{p_value}'

                result_file.write("\n" + str(csv_row))

        result_file.close()


DIGITS_OF_PRECISION = 4


def round_result(result):
    """ Round result by customized number of digits of precision according to the value of the result """
    return result
    # try:
    #     if (type(result) == np.float32 or type(result) == np.float64) and result != 0.0:
    #         if result >= 1.0:
    #             digits_of_precision = DIGITS_OF_PRECISION
    #         else:
    #             digits_of_precision = int(-np.log10(np.abs(result))) + DIGITS_OF_PRECISION
    #
    #         return np.round(result, digits_of_precision)
    #     else:
    #         return result
    # except:
    #     return result
