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
from SplinkStates.client import SplinkClient
from SplinkStates.server import SplinkServer
from Utils.utils import share_attrs, load_attrs


@app_state('Allele_Name', Role.BOTH)
class AlleleName(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Aggregate_Allele_Name', Role.COORDINATOR)
        self.register_transition('Non_Missing_counts', Role.PARTICIPANT)

    def run(self) -> str or None:
        """
        Sharing:
        allele_names
        sample_count
        Returns
        -------

        """
        snp_id_values = self.await_data()
        load_attrs(self)
        self.log(f"{len(snp_id_values)} SNPs are common among all clients")
        allele_names = self.allele_name_step(snp_id_values)
        self.send_data_to_coordinator(data=allele_names, use_smpc=False)
        share_attrs(self)
        if self.is_coordinator:
            self.store('allele_names_clients', self.gather_data())
        else:
            sleep(5)
        # self.send_data_to_coordinator(data=self.load('sample_count'), use_smpc=self.load('smpc_used'))
        self.send_data_to_coordinator(data=self.sample_count, use_smpc=self.load('smpc_used'))
        if self.is_coordinator:
            return 'Aggregate_Allele_Name'
        return 'Non_Missing_counts'

    def allele_name_step(self, snp_id_values):
        """ Initialize SNP (ID) values and first/second allele names based on global (common) SNP IDs first,
            and then, share allele names with the server """
        # update snp_id_values, snp_values, first/second allele names
        # based on the global SNP IDs by excluding those that do not exist in the other clients
        # The above-mentioned attributed are converted to list

        # update SNP values
        snp_values = list()
        # self.log(set(self.snp_values.keys()))
        for snp_id in snp_id_values:
            # snp_values.append(self.load('snp_values')[snp_id])
            snp_values.append(self.snp_values[snp_id])
        self.snp_values = snp_values

        # update first/second allele names and initialize allele names (shared with server)
        first_allele_names, second_allele_names, allele_names = [], [], [[], []]

        for snp_id in snp_id_values:
            first = self.first_allele_names[snp_id]
            second = self.second_allele_names[snp_id]
            first_allele_names.append(first)
            second_allele_names.append(second)

            # sort allele names to prevent revealing which allele is minor or major to the server
            if first < second:
                allele_names[0].append(first)
                allele_names[1].append(second)
            else:
                allele_names[0].append(second)
                allele_names[1].append(first)

        # self.store('first_allele_names', first_allele_names)
        # self.store('second_allele_names', second_allele_names)

        self.first_allele_names = first_allele_names
        self.second_allele_names = second_allele_names
        return allele_names


@app_state('Aggregate_Allele_Name', Role.COORDINATOR)
class AggregateAlleleName(AppState, SplinkServer):
    def __init__(self):
        SplinkServer.__init__(self)

    def register(self):
        self.register_transition('Non_Missing_counts', Role.COORDINATOR)

    def run(self) -> str or None:
        """
        Aggregating:
        allele_counts
        sample_count
        Returns
        -------

        """
        load_attrs(self)
        self.aggregate_alleles(self.load('allele_names_clients'))
        sample_count = self.aggregate_data(operation=SMPCOperation.ADD, use_smpc=self.load('smpc_used'))

        data_to_send = self.setup_next_chunk()
        self.broadcast_data(data_to_send)
        # self.broadcast_data(data=[self.current_chunk,
        #                           self.load('total_chunks'),
        #                           self.considered_snp_indices,
        #                           self.chunk_start_index,
        #                           self.chunk_end_index])
        share_attrs(self)
        return 'Non_Missing_counts'

    def aggregate_alleles(self, allele_names_clients):
        # initialize allele names for each SNP
        # for snp_index in np.arange(len(self.load('snp_id_values'))):
        for snp_index in np.arange(len(self.snp_id_values)):
            snp_alleles = [name[i][snp_index] for name in allele_names_clients for i in range(2)]
            self.allele_names[snp_index] = np.sort(np.unique(snp_alleles))

            # Ensure there are exactly two allele names for the SNP across all clients' datasets
            if len(self.allele_names[snp_index]) != 2:
                self.log(f"clients are using different allele names in their datasets!", LogLevel.FATAL)
                self.update(state=op_state.ERROR)
        self.attrs_to_share += ['allele_names']
