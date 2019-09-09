#
# ExpansionHunter Denovo
# Copyright 2016-2019 Illumina, Inc.
# All rights reserved.
#
# Author: Egor Dolzhenko <edolzhenko@illumina.com>,
#         Michael Eberle <meberle@illumina.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

import sys
import pytest

from core.wilcoxontest import *


class TestPermutationTest(object):
    def setup_method(self, test_method):
        self.cases = [8.50, 9.48, 8.65, 8.16, 8.83, 7.76, 8.63]
        self.controls = [8.27, 8.20, 8.25, 8.14, 9.00, 8.10, 7.20, 8.32, 7.70]
        self.expected_case_rank_sum = 75
        self.expected_pvalue = 0.057

    def test_permutation_pvalue(self):
        pvalue = calculate_permutation_pvalue(
            self.cases, self.controls, num_permutations=100000
        )

        assert abs(pvalue - self.expected_pvalue) < 0.01


class TestPermutationTestWithTies(object):
    def setup_method(self, test_method):
        self.cases = [0.45, 0.50, 0.61, 0.63, 0.75, 0.85, 0.93]
        self.controls = [0.44, 0.45, 0.52, 0.53, 0.56, 0.58, 0.58, 0.65, 0.79]
        self.expected_case_rank_sum = 71.5
        self.expected_pvalue = 0.105

    def test_permutation_pvalue(self):
        pvalue = calculate_permutation_pvalue(
            self.cases, self.controls, num_permutations=100000
        )

        assert abs(pvalue - self.expected_pvalue) < 0.01

    def test_approximate_pvalue(self):
        pvalue = calculate_approximate_pvalue(self.cases, self.controls)
        assert abs(pvalue - self.expected_pvalue) < 0.01
