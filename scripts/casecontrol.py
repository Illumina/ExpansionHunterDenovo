#!/usr/bin/env python3
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

import argparse
import numpy as np

import core.common as common
import casecontrol.locusworkflow
import casecontrol.motifworkflow


def initialize_parser():
    parser = argparse.ArgumentParser(
        prog="casecontrol", description="Case-control analysis of STR profiles"
    )
    return parser


def initialize_subparsers(parser):
    subparsers = parser.add_subparsers(help="command help")
    subparsers.required = True
    subparsers.dest = "command"
    return subparsers


def add_locus_command(subparsers):
    command_parser = subparsers.add_parser(
        "locus", help="Perform locus-based case-control analysis"
    )

    # Add required arguments
    required_args = command_parser.add_argument_group("required arguments")
    help = "TSV file describing all STR profiles"
    required_args.add_argument("--manifest", help=help, required=True)

    help = "JSON file with combined counts of anchored in-repeat reads"
    required_args.add_argument("--multisample-profile", help=help, required=True)

    help = "TSV file with results of the case-control analysis"
    required_args.add_argument("--output", help=help, required=True)

    # Add optional arguments
    min_irr_count = 5
    help = "Require regions to have at least this many in-repeat reads (default: {})"
    help = help.format(min_irr_count)
    command_parser.add_argument(
        "--min-inrepeat-reads", help=help, default=min_irr_count, type=int
    )

    help = "BED file with regions to which analysis should be restricted"
    command_parser.add_argument("--target-regions", help=help, default=None)

    wilcoxon_method = "normal"
    help = "Method of calculating Wilcoxon Rank-Sum Test p-value (default: {})"
    help = help.format(wilcoxon_method)
    command_parser.add_argument("--test-params", help=help, default=wilcoxon_method)

    # TODO: Implement logic for "resample(10000)" parameter
    # num_resamples = 1000000
    # help = "Number of iterations for the resampling test (default: {})"
    # help = help.format(num_resamples)

    return command_parser


def add_motif_command(subparsers):
    command_parser = subparsers.add_parser(
        "motif", help="Perform motif-based case-control analysis"
    )

    # Add required arguments
    required_args = command_parser.add_argument_group("required arguments")
    help = "TSV file describing all STR profiles"
    required_args.add_argument("--manifest", help=help, required=True)

    help = "JSON file with combined counts of anchored in-repeat reads"
    required_args.add_argument("--multisample-profile", help=help, required=True)

    help = "TSV file with results of the case-control analysis"
    required_args.add_argument("--output", help=help, required=True)

    # Add optional arguments
    min_irr_pair_count = 5
    help = "Require at least this many in-repeat read pairs for a given motif (default: {})"
    help = help.format(min_irr_pair_count)
    command_parser.add_argument(
        "--min-inrepeat-read-pairs", help=help, default=min_irr_pair_count, type=int
    )

    wilcoxon_method = "normal"
    help = "Method of calculating Wilcoxon Rank-Sum Test p-value (default: {})"
    help = help.format(wilcoxon_method)
    command_parser.add_argument("--test-params", help=help, default=wilcoxon_method)

    return command_parser


def decode_test_params(params):
    if params == "normal":
        return ("normal",)
    elif "permute_" in params:
        num_perms = int(params.replace("permute_", ""))
        return ("permute", num_perms)
    else:
        raise Exception("Unknown test parameters: {}".format(params))


def run_locus_workflow(args):
    params = casecontrol.locusworkflow.Parameters(
        manifest_path=args.manifest,
        multisample_profile_path=args.multisample_profile,
        min_inrepeat_reads=args.min_inrepeat_reads,
        output_path=args.output,
        target_region_path=args.target_regions,
        test_params=decode_test_params(args.test_params),
    )

    casecontrol.locusworkflow.run(params)


def run_motif_workflow(args):
    params = casecontrol.motifworkflow.Parameters(
        manifest_path=args.manifest,
        multisample_profile_path=args.multisample_profile,
        min_inrepeat_read_pairs=args.min_inrepeat_read_pairs,
        output_path=args.output,
        test_params=decode_test_params(args.test_params),
    )

    casecontrol.motifworkflow.run(params)


def main():
    np.random.seed(42)
    common.init_logger()
    parser = initialize_parser()
    subparsers = initialize_subparsers(parser)
    locus_command_parser = add_locus_command(subparsers)
    locus_command_parser.set_defaults(run_workflow=run_locus_workflow)
    motif_command_parser = add_motif_command(subparsers)
    motif_command_parser.set_defaults(run_workflow=run_motif_workflow)

    args = parser.parse_args()
    args.run_workflow(args)


if __name__ == "__main__":
    main()
