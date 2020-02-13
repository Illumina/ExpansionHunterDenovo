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

import json
import logging
from collections import namedtuple

import core.common as common

Parameters = namedtuple(
    "Parameters", ["manifest_path", "multisample_profile_path", "output_path"]
)


def run(params):
    with open(params.multisample_profile_path, "r") as profile_file:
        multisample_profile = json.load(profile_file)
    count_table = common.generate_table_with_irr_pair_counts(
        multisample_profile["Counts"]
    )

    logging.info("Loaded %i regions", len(count_table))

    logging.info("Normalizing counts")
    sample_stats = multisample_profile["Parameters"]
    common.depth_normalize_counts(sample_stats, count_table)

    manifest = common.load_manifest(params.manifest_path)
    sample_status = common.extract_case_control_assignments(manifest)

    header = "motif\ttop_case_zscore\thigh_case_counts\tcounts"
    with open(params.output_path, "w") as results_file:
        print(header, file=results_file)
        for row in count_table:
            top_case_zscore, cases_with_high_counts = common.run_zscore_analysis(
                sample_status, row["sample_counts"]
            )

            if len(cases_with_high_counts) == 0:
                continue

            encoded_case_info = ",".join(
                "{}:{:.2f}".format(s, c) for s, c in cases_with_high_counts.items()
            )
            count_encoding = ",".join(
                ["{:.2f}".format(c) for _, c in row["sample_counts"].items()]
            )

            motif = row["unit"]
            print(
                motif,
                "{:.2f}".format(top_case_zscore),
                encoded_case_info,
                count_encoding,
                sep="\t",
                file=results_file,
            )
