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

from FeatureCloud.app.engine.app import app_state, AppState, Role, LogLevel
from FeatureCloud.app.engine.app import State as op_state
from CustomStates import ConfigState
from Utils.gwas_dataset import GwasDataset
import numpy as np
from States import ALGORITHM
from SplinkStates.client import SplinkClient
from SplinkStates.server import SplinkServer
from Utils.utils import share_attrs, load_attrs

name = 'splink'


@app_state(name='initial', role=Role.BOTH, app_name=name)
class Init(ConfigState.State, SplinkClient):
    def __init__(self, app_name, input_dir: str = "/mnt/input", output_dir: str = "/mnt/output"):
        ConfigState.State.__init__(self, app_name, input_dir, output_dir)
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_SNP', Role.COORDINATOR)
        self.register_transition('Allele_Name', Role.PARTICIPANT)

    def run(self) -> str or None:
        self.lazy_init()
        self.read_config()
        self.finalize_config()
        self.config['covariates'] = () if self.config['covariates'] is None else \
            tuple(filter(None, self.config['covariates'].split(',')))
        self.config['chunk_size'] = self.config['chunk_size'] * 1000
        self.store('config', self.config)
        self.store('smpc_used', self.config.get('use_smpc', False))
        self.read_dataset()
        non_zero_snp_ids = self.snp_id_step()
        self.send_data_to_coordinator(data=non_zero_snp_ids, use_smpc=False)

        # self.store_attrs()
        share_attrs(self)

        if self.is_coordinator:
            return 'Aggregate_SNP'
        return 'Allele_Name'
        # return 'Allele_Name'


    def read_dataset(self):
        gwas_dataset = GwasDataset(bed_file_path=self.load('input_files')['data'][0],
                                   phenotype_file_path=self.load('input_files')['phenotype'][0],
                                   covariate_file_path=self.load('input_files')['covariate'][0],
                                   phenotype_name=self.config['phenotype_name'],
                                   covariate_names=self.config['covariate_names'])

        self.log("Opening and pre-processing the GWAS dataset ...")
        gwas_dataset.open_and_preprocess()
        # if gwas_dataset.is_operation_failed():
        #     self.log(gwas_dataset.get_error_message(), LogLevel.FATAL)
        #     self.update(state=op_state.ERROR)

        # phenotypes should be binary for logistic regression and chi-square
        if not gwas_dataset.is_phenotype_binary() and \
                (self.load('config')['algorithm'] in [ALGORITHM.CHI_SQUARE, ALGORITHM.LINEAR_REGRESSION]):
            self.log(f"Phenotype values must be binary for {self.load('config')['algorithm']} tests!",
                     LogLevel.FATAL)
            self.update(state=op_state.ERROR)

        # phenotype values should be quantitative for linear regression
        if gwas_dataset.is_phenotype_binary() and self.load('config')['algorithm'] == ALGORITHM.LINEAR_REGRESSION:
            self.log(f"Phenotype values must be quantitative for {self.load('config')['algorithm']} tests!",
                     LogLevel.FATAL)
            self.update(state=op_state.ERROR)

        # log general info about the gwas dataset
        self.log(gwas_dataset.get_dataset_info())

        # #### initialize attributes required to compute local parameters

        # initialize fam file related attributes
        self.sex_values = gwas_dataset.get_sex_values()
        self.phenotype_values = gwas_dataset.get_phenotype_values()
        self.sample_count = gwas_dataset.get_sample_count()

        # initialize bim file related attributes
        self.snp_id_values = gwas_dataset.get_snp_id_values()
        self.first_allele_names = gwas_dataset.get_first_allele_names()
        self.second_allele_names = gwas_dataset.get_second_allele_names()

        # initialize bed file related attribute
        self.snp_values = gwas_dataset.get_snp_values()

        # initialize covariate file related attributes
        self.covariate_values = gwas_dataset.get_covariate_values()

        # initialize attributes to speed-up the creation of the feature matrix
        self.non_missing_index_values = gwas_dataset.get_non_missing_index_values()
        self.covariate_matrix = gwas_dataset.get_covariate_matrix()

    def snp_id_step(self):
        """ share SNP IDs with the server """

        # share the SNP IDs whose minor allele frequency is non-zero with the server
        non_zero_snp_ids = np.array([snp_id for snp_id in self.snp_id_values if self.first_allele_names[snp_id] != '0'])
        # non_zero_snp_ids = [snp_id for snp_id in self.snp_id_values if self.load('first_allele_names')[snp_id] != '0']
        return non_zero_snp_ids
        # return non_zero_snp_ids.tolist()


@app_state('Aggregate_SNP', Role.COORDINATOR)
class AggregateSNP(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Allele_Name', Role.COORDINATOR)

    def run(self) -> str or None:
        snp_ids_clients = self.await_data(unwrap=False)
        self.snp_id_step(snp_ids_clients)
        self.broadcast_data(self.snp_id_values)
        share_attrs(self)
        return 'Allele_Name'

    def snp_id_step(self, snp_ids_clients):
        """ Compute the intersection of the SNP IDs from the clients """

        # compute the SNP IDs common among all clients
        self.log("Inside SNP ID Agg")
        intersect_snp_ids = set(snp_ids_clients[0])
        for client_snp_ids in snp_ids_clients:
            intersect_snp_ids = intersect_snp_ids.intersection(client_snp_ids)
        self.snp_id_values = np.array(list(intersect_snp_ids), dtype="S")
        if len(self.snp_id_values) == 0:
            self.log("There is no SNP common among all clients!", LogLevel.FATAL)
            self.update(state=op_state.ERROR)

        # initialize chunks
        self.init_chunks()

    def init_chunks(self):
        """ Set the total number of chunks and start/end indices of the chunks """
        self.total_chunks = int(np.ceil(len(self.snp_id_values) / self.load('config')['chunk_size']))
        for split in np.array_split(np.arange(len(self.snp_id_values)), self.total_chunks):
            self.start_indices_chunks.append(split[0])
            self.end_indices_chunks.append(split[-1] + 1)
        self.attrs_to_share += ['start_indices_chunks', 'end_indices_chunks']
