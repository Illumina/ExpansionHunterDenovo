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

import numpy as np
import scipy.special
import scipy.stats


def calculate_approximate_pvalue(cases, controls):
    all_counts = cases + controls
    ranks = scipy.stats.rankdata(all_counts)
    case_rank_sum = np.sum(ranks[: len(cases)])
    num_counts = len(cases) + len(controls)

    mu_cases = len(cases) * (num_counts + 1) / 2
    sigma_cases = np.sqrt(len(cases) * len(controls) * (num_counts + 1) / 12)
    z_cases = (case_rank_sum - mu_cases) / sigma_cases

    return 1 - scipy.stats.norm.cdf(z_cases)


def calculate_permutation_pvalue(cases, controls, num_permutations):
    all_counts = cases + controls
    ranks = scipy.stats.rankdata(all_counts)
    num_cases = len(cases)
    true_case_rank_sum = np.sum(ranks[:num_cases])

    permuted_case_ranks = np.random.choice(ranks, size=(num_permutations, num_cases))
    permuted_case_rank_sums = np.sum(permuted_case_ranks, axis=1)

    num_case_rank_sums_as_extreme_as_true = np.sum(
        permuted_case_rank_sums >= true_case_rank_sum
    )

    return (num_case_rank_sums_as_extreme_as_true + 1) / (num_permutations + 1)


def wilcoxon_rank_sum_test(test_params, cases, controls):
    method, *params = test_params
    if method == "normal":
        return calculate_approximate_pvalue(cases, controls)
    elif method == "permute":
        num_perms = params[0]
        return calculate_permutation_pvalue(cases, controls, num_perms)
    else:
        assert False, "{} is an unknown method type".format(method)
