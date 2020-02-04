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

import collections
import logging
import json
import scipy.stats as stats
import numpy as np

from . import regiontools

from .wilcoxontest import wilcoxon_rank_sum_test


def init_logger():
    logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.INFO)


def load_manifest(manifest_path):
    """Extract sample information from a manifest file.

    """

    # pylint: disable=I0011,C0103
    Sample = collections.namedtuple("Sample", "id status path")

    samples = []
    with open(manifest_path, "r") as manifest_file:
        for line in manifest_file:
            sample_id, status, path = line.split()

            if status not in ["case", "control"]:
                message = (
                    'Sample status must be either "case" or "control";'
                    ' instead got "{}"'
                )
                raise Exception(message.format(status))

            sample = Sample(id=sample_id, status=status, path=path)
            samples.append(sample)
    return samples


def filter_counts_by_magnitude(count_table, count_cutoff):
    filtered_count_table = []
    for row in count_table:
        max_count = max(count for _, count in row["sample_counts"].items())
        if max_count >= count_cutoff:
            filtered_count_table.append(row)

    return filtered_count_table


def filter_counts_by_region(count_table, target_regions):
    filtered_count_table = []
    for row in count_table:
        region_encoding = row["region"]
        chrom, coords = region_encoding.rsplit(":", 1)
        start, end = coords.split("-")
        start, end = int(start), int(end)
        region = regiontools.Region(chrom, start, end)
        overlaps_target_region = any(
            regiontools.compute_distance(region, target) == 0
            for target in target_regions
        )

        if overlaps_target_region:
            filtered_count_table.append(row)

    return filtered_count_table


def extract_case_control_assignments(samples):
    sample_status = {}
    for sample in samples:
        sample_status[sample.id] = sample.status
    return sample_status


def test_samples(test_params, sample_status, sample_counts):
    control_samples = [
        sample for sample, status in sample_status.items() if status == "control"
    ]
    case_samples = [
        sample for sample, status in sample_status.items() if status == "case"
    ]

    control_counts = [
        sample_counts[s] if s in sample_counts else 0 for s in control_samples
    ]

    case_counts = [sample_counts[s] if s in sample_counts else 0 for s in case_samples]

    pvalue = wilcoxon_rank_sum_test(test_params, case_counts, control_counts)

    return pvalue


def compare_counts(test_params, sample_status, count_table):
    for row in count_table:
        # Generate counts before testing
        pvalue = test_samples(test_params, sample_status, row["sample_counts"])
        row["pvalue"] = pvalue


def correct_pvalues(count_table):
    num_tests = len(count_table)
    for row in count_table:
        row["bonf_pvalue"] = min(row["pvalue"] * num_tests, 1.0)


def normalize_count(sample_depth, count, target_depth=40):
    return target_depth * count / sample_depth


def depth_normalize_counts(sample_stats, count_table):
    depths = sample_stats["Depths"]
    for row in count_table:
        row["sample_counts"] = {
            s: normalize_count(depths[s], c) for s, c in row["sample_counts"].items()
        }


def generate_table_with_irr_pair_counts(combined_counts):
    count_table = []
    for unit, rec in combined_counts.items():
        if "IrrPairCounts" not in rec:
            continue

        sample_counts = rec["IrrPairCounts"]
        table_row = {"unit": unit, "sample_counts": sample_counts}
        count_table.append(table_row)

    return count_table


def resample_quantiles(counts, num_resamples, target_quantile_value):
    resamples = np.random.choice(counts, len(counts) * num_resamples)
    resamples = np.split(resamples, num_resamples)

    resampled_quantiles = []
    for resample in resamples:
        quantile = np.quantile(resample, target_quantile_value)
        resampled_quantiles.append(quantile)

    return resampled_quantiles


def run_zscore_analysis(sample_status, sample_counts):
    raw_counts = [sample_counts.get(sample, 0) for sample, _ in sample_status.items()]
    quantiles = resample_quantiles(raw_counts, 100, 0.95)
    (mu, sigma) = stats.norm.fit(quantiles)
    sigma = max(sigma, 1)

    case_counts = {
        sample: sample_counts.get(sample, 0)
        for sample, status in sample_status.items()
        if status == "case"
    }

    assert len(case_counts) >= 1, "Manifest must contain at least one case"

    cases_with_high_counts = {}
    top_zscore = -1
    for sample, count in case_counts.items():
        zscore = (count - mu) / sigma

        if zscore > 1.0:
            cases_with_high_counts[sample] = count
            top_zscore = max(top_zscore, zscore)

    return (top_zscore, cases_with_high_counts)
