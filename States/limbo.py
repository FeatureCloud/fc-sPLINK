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


@app_state('Limbo', Role.PARTICIPANT)
class Limbo(AppState, SplinkClient):
    def __init__(self):
        SplinkClient.__init__(self)

    def register(self):
        self.register_transition('Non_Missing_counts', Role.PARTICIPANT)
        self.register_transition('terminal', Role.PARTICIPANT)

    def run(self) -> str or None:
        data = self.await_data()
        load_attrs(self)
        if data is None:
            return 'terminal'
        self.current_chunk, self.total_chunks, self.snp_indices, self.chunk_start_index, self.chunk_end_index = data
        share_attrs(self)
        return 'Non_Missing_counts'
