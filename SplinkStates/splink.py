class Splink:
    def __init__(self):
        self._allele_names = dict()
        self._non_missing_sample_counts = {}
        self._allele_counts = {}
        self._minor_allele_names = {}
        self._major_allele_names = {}
        self._minor_allele_counts = {}
        self._major_allele_counts = {}
        self._minor_allele_frequencies = {}
        self._major_allele_frequencies = {}
        self._contingency_tables = {}
        self._maf_case = {}
        self._maf_control = {}
        self._chi_square_values = {}
        self._odd_ratio_values = {}
        self._xt_x_matrices = {}
        self._xt_y_vectors = {}
        self._xt_x_inverse_matrices = {}
        self._sse_values = {}
        self._gradient_vectors = {}
        self._hessian_matrices = {}
        self._new_log_likelihood_values = {}
        self._new_beta_values = {}
        self._log_likelihood_values = {}
        self._beta_values = {}
        self._std_error_values = {}
        self._t_stat_values = {}
        self._p_values = {}
        self._sex_values = []
        self._phenotype_values = []
        self._sample_count = 0
        self._snp_id_values = []
        self._first_allele_names = {}
        self._second_allele_names = {}
        self._snp_values = {}
        self._covariate_values = {}
        self._non_missing_index_values = []
        self._covariate_matrix = []
        self._current_chunk = 0
        self._total_chunks = -1
        self._sub_chunk_start_indices = []
        self._sub_chunk_end_indices = []
        self._snp_indices = set()
        self._current_chunk_size = 0
        self._xt_x_matrix = {}
        self._xt_y_vector = {}
        self._current_beta_iteration = -1
        self._chromosome_number_all_chunks = []
        self._base_pair_distance_all_chunks = []
        self._p_value_all_chunks = []
        self._start_indices_chunks = []
        self._end_indices_chunks = []
        self._chunk_start_index = -1
        self._chunk_end_index = -1
        self._considered_snp_indices = set()
        self._in_process_snp_indices = set()
        self._considered_in_process_snp_indices = set()
        self.attrs_to_share = []
        self.delta_log_likelihood_threshold = 0.0001


    @property
    def allele_names(self):
        return self._allele_names

    @allele_names.setter
    def allele_names(self, value):
        self._allele_names = value
        self.attrs_to_share.append('allele_names')

    @property
    def non_missing_sample_counts(self):
        return self._non_missing_sample_counts

    @non_missing_sample_counts.setter
    def non_missing_sample_counts(self, value):
        self._non_missing_sample_counts = value
        self.attrs_to_share.append('non_missing_sample_counts')

    @property
    def allele_counts(self):
        return self._allele_counts

    @allele_counts.setter
    def allele_counts(self, value):
        self._allele_counts = value
        self.attrs_to_share.append('allele_counts')

    @property
    def minor_allele_names(self):
        return self._minor_allele_names

    @minor_allele_names.setter
    def minor_allele_names(self, value):
        self._minor_allele_names = value
        self.attrs_to_share.append('minor_allele_names')

    @property
    def major_allele_names(self):
        return self._major_allele_names

    @major_allele_names.setter
    def major_allele_names(self, value):
        self._major_allele_names = value
        self.attrs_to_share.append('major_allele_names')

    @property
    def minor_allele_counts(self):
        return self._minor_allele_counts

    @minor_allele_counts.setter
    def minor_allele_counts(self, value):
        self._minor_allele_counts = value
        self.attrs_to_share.append('minor_allele_counts')

    @property
    def major_allele_counts(self):
        return self._major_allele_counts

    @major_allele_counts.setter
    def major_allele_counts(self, value):
        self._major_allele_counts = value
        self.attrs_to_share.append('major_allele_counts')

    @property
    def minor_allele_frequencies(self):
        return self._minor_allele_frequencies

    @minor_allele_frequencies.setter
    def minor_allele_frequencies(self, value):
        self._minor_allele_frequencies = value
        self.attrs_to_share.append('minor_allele_frequencies')

    @property
    def major_allele_frequencies(self):
        return self._major_allele_frequencies

    @major_allele_frequencies.setter
    def major_allele_frequencies(self, value):
        self._major_allele_frequencies = value
        self.attrs_to_share.append('major_allele_frequencies')

    @property
    def contingency_tables(self):
        return self._contingency_tables

    @contingency_tables.setter
    def contingency_tables(self, value):
        self._contingency_tables = value
        self.attrs_to_share.append('contingency_tables')

    @property
    def maf_case(self):
        return self._maf_case

    @maf_case.setter
    def maf_case(self, value):
        self._maf_case = value
        self.attrs_to_share.append('maf_case')

    @property
    def maf_control(self):
        return self._maf_control

    @maf_control.setter
    def maf_control(self, value):
        self._maf_control = value
        self.attrs_to_share.append('maf_control')

    @property
    def chi_square_values(self):
        return self._chi_square_values

    @chi_square_values.setter
    def chi_square_values(self, value):
        self._chi_square_values = value
        self.attrs_to_share.append('chi_square_values')

    @property
    def odd_ratio_values(self):
        return self._odd_ratio_values

    @odd_ratio_values.setter
    def odd_ratio_values(self, value):
        self._odd_ratio_values = value
        self.attrs_to_share.append('odd_ratio_values')

    @property
    def xt_x_matrices(self):
        return self._xt_x_matrices

    @xt_x_matrices.setter
    def xt_x_matrices(self, value):
        self._xt_x_matrices = value
        self.attrs_to_share.append('xt_x_matrices')

    @property
    def xt_y_vectors(self):
        return self._xt_y_vectors

    @xt_y_vectors.setter
    def xt_y_vectors(self, value):
        self._xt_y_vectors = value
        self.attrs_to_share.append('xt_y_vectors')

    @property
    def xt_x_inverse_matrices(self):
        return self._xt_x_inverse_matrices

    @xt_x_inverse_matrices.setter
    def xt_x_inverse_matrices(self, value):
        self._xt_x_inverse_matrices = value
        self.attrs_to_share.append('xt_x_inverse_matrices')

    @property
    def sse_values(self):
        return self._sse_values

    @sse_values.setter
    def sse_values(self, value):
        self._sse_values = value
        self.attrs_to_share.append('sse_values')

    @property
    def gradient_vectors(self):
        return self._gradient_vectors

    @gradient_vectors.setter
    def gradient_vectors(self, value):
        self._gradient_vectors = value
        self.attrs_to_share.append('gradient_vectors')

    @property
    def hessian_matrices(self):
        return self._hessian_matrices

    @hessian_matrices.setter
    def hessian_matrices(self, value):
        self._hessian_matrices = value
        self.attrs_to_share.append('hessian_matrices')

    @property
    def new_log_likelihood_values(self):
        return self._new_log_likelihood_values

    @new_log_likelihood_values.setter
    def new_log_likelihood_values(self, value):
        self._new_log_likelihood_values = value
        self.attrs_to_share.append('new_log_likelihood_values')

    @property
    def new_beta_values(self):
        return self._new_beta_values

    @new_beta_values.setter
    def new_beta_values(self, value):
        self._new_beta_values = value
        self.attrs_to_share.append('new_beta_values')

    @property
    def log_likelihood_values(self):
        return self._log_likelihood_values

    @log_likelihood_values.setter
    def log_likelihood_values(self, value):
        self._log_likelihood_values = value
        self.attrs_to_share.append('log_likelihood_values')

    @property
    def beta_values(self):
        return self._beta_values

    @beta_values.setter
    def beta_values(self, value):
        self._beta_values = value
        self.attrs_to_share.append('beta_values')

    @property
    def std_error_values(self):
        return self._std_error_values

    @std_error_values.setter
    def std_error_values(self, value):
        self._std_error_values = value
        self.attrs_to_share.append('std_error_values')

    @property
    def t_stat_values(self):
        return self._t_stat_values

    @t_stat_values.setter
    def t_stat_values(self, value):
        self._t_stat_values = value
        self.attrs_to_share.append('t_stat_values')

    @property
    def p_values(self):
        return self._p_values

    @p_values.setter
    def p_values(self, value):
        self._p_values = value
        self.attrs_to_share.append('p_values')

    @property
    def sex_values(self):
        return self._sex_values

    @sex_values.setter
    def sex_values(self, value):
        self._sex_values = value
        self.attrs_to_share.append('sex_values')

    @property
    def phenotype_values(self):
        return self._phenotype_values

    @phenotype_values.setter
    def phenotype_values(self, value):
        self._phenotype_values = value
        self.attrs_to_share.append('phenotype_values')

    @property
    def sample_count(self):
        return self._sample_count

    @sample_count.setter
    def sample_count(self, value):
        self._sample_count = value
        self.attrs_to_share.append('sample_count')

    @property
    def snp_id_values(self):
        return self._snp_id_values

    @snp_id_values.setter
    def snp_id_values(self, value):
        self._snp_id_values = value
        self.attrs_to_share.append('snp_id_values')

    @property
    def first_allele_names(self):
        return self._first_allele_names

    @first_allele_names.setter
    def first_allele_names(self, value):
        self._first_allele_names = value
        self.attrs_to_share.append('first_allele_names')

    @property
    def second_allele_names(self):
        return self._second_allele_names

    @second_allele_names.setter
    def second_allele_names(self, value):
        self._second_allele_names = value
        self.attrs_to_share.append('second_allele_names')

    @property
    def snp_values(self):
        return self._snp_values

    @snp_values.setter
    def snp_values(self, value):
        self._snp_values = value
        self.attrs_to_share.append('snp_values')

    @property
    def covariate_values(self):
        return self._covariate_values

    @covariate_values.setter
    def covariate_values(self, value):
        self._covariate_values = value
        self.attrs_to_share.append('covariate_values')

    @property
    def non_missing_index_values(self):
        return self._non_missing_index_values

    @non_missing_index_values.setter
    def non_missing_index_values(self, value):
        self._non_missing_index_values = value
        self.attrs_to_share.append('non_missing_index_values')

    @property
    def covariate_matrix(self):
        return self._covariate_matrix

    @covariate_matrix.setter
    def covariate_matrix(self, value):
        self._covariate_matrix = value
        self.attrs_to_share.append('covariate_matrix')

    @property
    def current_chunk(self):
        return self._current_chunk

    @current_chunk.setter
    def current_chunk(self, value):
        self._current_chunk = value
        self.attrs_to_share.append('current_chunk')

    @property
    def total_chunks(self):
        return self._total_chunks

    @total_chunks.setter
    def total_chunks(self, value):
        self._total_chunks = value
        self.attrs_to_share.append('total_chunks')

    @property
    def sub_chunk_start_indices(self):
        return self._sub_chunk_start_indices

    @sub_chunk_start_indices.setter
    def sub_chunk_start_indices(self, value):
        self._sub_chunk_start_indices = value
        self.attrs_to_share.append('sub_chunk_start_indices')

    @property
    def sub_chunk_end_indices(self):
        return self._sub_chunk_end_indices

    @sub_chunk_end_indices.setter
    def sub_chunk_end_indices(self, value):
        self._sub_chunk_end_indices = value
        self.attrs_to_share.append('sub_chunk_end_indices')

    @property
    def snp_indices(self):
        return self._snp_indices

    @snp_indices.setter
    def snp_indices(self, value):
        self._snp_indices = value
        self.attrs_to_share.append('snp_indices')

    @property
    def current_chunk_size(self):
        return self._current_chunk_size

    @current_chunk_size.setter
    def current_chunk_size(self, value):
        self._current_chunk_size = value
        self.attrs_to_share.append('current_chunk_size')

    @property
    def xt_x_matrix(self):
        return self._xt_x_matrix

    @xt_x_matrix.setter
    def xt_x_matrix(self, value):
        self._xt_x_matrix = value
        self.attrs_to_share.append('xt_x_matrix')

    @property
    def xt_y_vector(self):
        return self._xt_y_vector

    @xt_y_vector.setter
    def xt_y_vector(self, value):
        self._xt_y_vector = value
        self.attrs_to_share.append('xt_y_vector')

    @property
    def current_beta_iteration(self):
        return self._current_beta_iteration

    @current_beta_iteration.setter
    def current_beta_iteration(self, value):
        self._current_beta_iteration = value
        self.attrs_to_share.append('current_beta_iteration')

    @property
    def chromosome_number_all_chunks(self):
        return self._chromosome_number_all_chunks

    @chromosome_number_all_chunks.setter
    def chromosome_number_all_chunks(self, value):
        self._chromosome_number_all_chunks = value
        self.attrs_to_share.append('chromosome_number_all_chunks')

    @property
    def base_pair_distance_all_chunks(self):
        return self._base_pair_distance_all_chunks

    @base_pair_distance_all_chunks.setter
    def base_pair_distance_all_chunks(self, value):
        self._base_pair_distance_all_chunks = value
        self.attrs_to_share.append('base_pair_distance_all_chunks')

    @property
    def p_value_all_chunks(self):
        return self._p_value_all_chunks

    @p_value_all_chunks.setter
    def p_value_all_chunks(self, value):
        self._p_value_all_chunks = value
        self.attrs_to_share.append('p_value_all_chunks')

    @property
    def start_indices_chunks(self):
        return self._start_indices_chunks

    @start_indices_chunks.setter
    def start_indices_chunks(self, value):
        self._start_indices_chunks = value
        self.attrs_to_share.append('start_indices_chunks')

    @property
    def end_indices_chunks(self):
        return self._end_indices_chunks

    @end_indices_chunks.setter
    def end_indices_chunks(self, value):
        self._end_indices_chunks = value
        self.attrs_to_share.append('end_indices_chunks')

    @property
    def chunk_start_index(self):
        return self._chunk_start_index

    @chunk_start_index.setter
    def chunk_start_index(self, value):
        self._chunk_start_index = value
        self.attrs_to_share.append('chunk_start_index')

    @property
    def chunk_end_index(self):
        return self._chunk_end_index

    @chunk_end_index.setter
    def chunk_end_index(self, value):
        self._chunk_end_index = value
        self.attrs_to_share.append('chunk_end_index')

    @property
    def considered_snp_indices(self):
        return self._considered_snp_indices

    @considered_snp_indices.setter
    def considered_snp_indices(self, value):
        self._considered_snp_indices = value
        self.attrs_to_share.append('considered_snp_indices')

    @property
    def in_process_snp_indices(self):
        return self._in_process_snp_indices

    @in_process_snp_indices.setter
    def in_process_snp_indices(self, value):
        self._in_process_snp_indices = value
        self.attrs_to_share.append('in_process_snp_indices')

    @property
    def considered_in_process_snp_indices(self):
        return self._considered_in_process_snp_indices

    @considered_in_process_snp_indices.setter
    def considered_in_process_snp_indices(self, value):
        self._considered_in_process_snp_indices = value
        self.attrs_to_share.append('considered_in_process_snp_indices')

