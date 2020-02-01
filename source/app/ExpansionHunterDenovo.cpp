//
// ExpansionHunter Denovo
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Michael Eberle <meberle@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include <iostream>

#include <boost/program_options.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "app/Version.hh"
#include "merge/MergeWorkflow.hh"
#include "profile/ProfileWorkflow.hh"

namespace po = boost::program_options;

using std::string;

enum class Workflow
{
    kProfile,
    kMerge
};

int readBaselineOptions(int argc, char** argv)
{
    // clang-format off
    string helpMessage = kProgramVersion
        + "\n\nUsage: ExpansionHunterDenovo <command> [options]\n\n"
        + "Commands:\n"
        + " profile  Compute genome-wide STR profile\n"
        + " merge    Generate multisample STR profile from single-sample profiles";

    po::options_description options("Available commands");
    options.add_options()
        ("help", "Print help message")
        ("version", "Print version");
    // clang-format on

    po::variables_map optionsMap;
    po::store(po::command_line_parser(argc, argv).options(options).run(), optionsMap);
    po::notify(optionsMap);

    if (optionsMap.count("help"))
    {
        std::cerr << helpMessage << std::endl;
        return 0;
    }

    if (optionsMap.count("version"))
    {
        std::cerr << kProgramVersion << std::endl;
        return 0;
    }

    std::cerr << helpMessage << std::endl;
    return 1;
}

int runProfileWorkflow(int argc, char** argv)
{
    string helpHeader = "Usage: ExpansionHunterDenovo profile [options]\n\n";

    string pathToReads;
    string pathToReference;
    string outputPrefix;
    int shortestUnitToConsider = 2;
    int longestUnitToConsider = 20;
    int minMapqOfAnchorRead = 50;
    int maxMapqOfInrepeatRead = 40;
    bool enableReadLog = false;

    // clang-format off
    po::options_description options("Available options");
    options.add_options()
        ("help", "Print help message")
        ("reads", po::value<string>(&pathToReads)->required(), "BAM or CRAM file with aligned reads")
        ("reference", po::value<string>(&pathToReference)->required(), "FASTA file with reference assembly")
        ("output-prefix", po::value<string>(&outputPrefix)->required(), "Prefix for the output files")
        ("min-unit-len", po::value<int>(&shortestUnitToConsider)->default_value(shortestUnitToConsider), "Shortest repeat unit to consider")
        ("max-unit-len", po::value<int>(&longestUnitToConsider)->default_value(longestUnitToConsider), "Longest repeat unit to consider")
        ("min-anchor-mapq", po::value<int>(&minMapqOfAnchorRead)->default_value(minMapqOfAnchorRead), "Minimum MAPQ of an anchor read")
        ("max-irr-mapq", po::value<int>(&maxMapqOfInrepeatRead)->default_value(maxMapqOfInrepeatRead), "Maximum MAPQ of an in-repeat read")
        ("log-reads", po::bool_switch(&enableReadLog), "Log informative reads");
    // clang-format on

    po::variables_map optionsMap;

    if (argc == 1)
    {
        std::cerr << helpHeader << options << std::endl;
        return 1;
    }

    try
    {
        po::store(po::command_line_parser(argc, argv).options(options).run(), optionsMap);

        if (optionsMap.count("help"))
        {
            std::cerr << helpHeader << options << std::endl;
            return 0;
        }

        po::notify(optionsMap);
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    spdlog::info("Starting {} profile workflow", kProgramVersion);

    Interval motifSizeRange(shortestUnitToConsider, longestUnitToConsider);
    ProfileWorkflowParameters params(
        outputPrefix, enableReadLog, pathToReads, pathToReference, motifSizeRange, minMapqOfAnchorRead,
        maxMapqOfInrepeatRead);

    return runProfileWorkflow(params);
}

int runMergeWorkflow(int argc, char** argv)
{
    string helpHeader = "Usage: ExpansionHunterDenovo merge [options]\n\n";

    string pathToReference;
    string pathToManifest;
    string outputPrefix;
    int shortestUnitToConsider = 2;
    int longestUnitToConsider = 20;

    // clang-format off
    po::options_description options("Available options");
    options.add_options()
        ("help", "Print help message")
        ("reference", po::value<string>(&pathToReference)->required(), "FASTA file with reference assembly")
        ("manifest", po::value<string>(&pathToManifest)->required(), "TSV with sample names and absolute paths")
        ("output-prefix", po::value<string>(&outputPrefix)->required(), "Prefix for the output files")
        ("min-unit-len", po::value<int>(&shortestUnitToConsider)->default_value(shortestUnitToConsider), "Shortest repeat unit to consider")
        ("max-unit-len", po::value<int>(&longestUnitToConsider)->default_value(longestUnitToConsider), "Longest repeat unit to consider");
    // clang-format on

    po::variables_map optionsMap;

    if (argc == 1)
    {
        std::cerr << helpHeader << options << std::endl;
        return 1;
    }

    try
    {
        po::store(po::command_line_parser(argc, argv).options(options).run(), optionsMap);

        if (optionsMap.count("help"))
        {
            std::cerr << helpHeader << options << std::endl;
            return 0;
        }

        po::notify(optionsMap);
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << std::endl;
        return 1;
    }

    spdlog::info("Starting {} profile workflow", kProgramVersion);

    MergeWorkflowParameters params(
        pathToReference, outputPrefix, pathToManifest, shortestUnitToConsider, longestUnitToConsider);
    return runMergeWorkflow(params);
}

int main(int argc, char** argv)
{
    try
    {
        if (argc == 1)
        {
            return readBaselineOptions(argc, argv);
        }

        const string command(argv[1]);
        if (command == "profile")
        {
            return runProfileWorkflow(argc - 1, argv + 1);
        }
        else if (command == "merge")
        {
            return runMergeWorkflow(argc - 1, argv + 1);
        }
        else
        {
            return readBaselineOptions(argc, argv);
        }
    }
    catch (const std::exception& exception)
    {
        spdlog::error(exception.what());
        return 1;
    }
}
