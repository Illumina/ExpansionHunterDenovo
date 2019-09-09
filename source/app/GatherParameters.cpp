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

#include "GatherParameters.hh"

#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>

#include "app/Version.hh"

using boost::optional;
using std::string;

namespace po = boost::program_options;

struct CommandLineParameters
{
    string pathToReads;
    string pathToReference;
    string outputPrefix;
    int shortestUnitToConsider = 2;
    int longestUnitToConsider = 20;
    int minMapqOfAnchorRead = 50;
    int maxMapqOfInrepeatRead = 40;
};

optional<CommandLineParameters> loadCommandLineParams(int argc, char** argv)
{
    CommandLineParameters params;
    // clang-format off
    po::options_description options("Available options");
    options.add_options()
        ("help", "Print help message")
        ("version", "Print version")
        ("reads", po::value<string>(&params.pathToReads)->required(), "BAM or CRAM file with aligned reads")
        ("reference", po::value<string>(&params.pathToReference)->required(), "FASTA file with reference assembly")
        ("output-prefix", po::value<string>(&params.outputPrefix)->required(), "Prefix for the output files")
        ("min-unit-len", po::value<int>(&params.shortestUnitToConsider)->default_value(2), "Shortest repeat unit to consider")
        ("max-unit-len", po::value<int>(&params.longestUnitToConsider)->default_value(20), "Longest repeat unit to consider")
        ("min-anchor-mapq", po::value<int>(&params.minMapqOfAnchorRead)->default_value(50), "Minimum MAPQ of an anchor read")
        ("max-irr-mapq", po::value<int>(&params.maxMapqOfInrepeatRead)->default_value(40), "Maximum MAPQ of an in-repeat read");
    // clang-format on

    po::variables_map optionsMap;
    po::store(po::command_line_parser(argc, argv).options(options).run(), optionsMap);

    if (argc == 1)
    {
        std::cerr << options << std::endl;
        throw std::runtime_error("No options provided");
    }

    if (optionsMap.count("help"))
    {
        std::cerr << options << std::endl;
        return boost::none;
    }

    if (optionsMap.count("version"))
    {
        std::cerr << kProgramVersion << std::endl;
        return boost::none;
    }

    po::notify(optionsMap);

    return params;
}

optional<ProgramParameters> loadParameters(int argc, char** argv)
{
    auto optionalCommandlineParams = loadCommandLineParams(argc, argv);

    if (optionalCommandlineParams)
    {
        PathParameters paths(
            optionalCommandlineParams->pathToReads, optionalCommandlineParams->pathToReference,
            optionalCommandlineParams->outputPrefix);

        const int readLength = 100;

        HeuristicParameters heuristics(
            optionalCommandlineParams->shortestUnitToConsider, optionalCommandlineParams->longestUnitToConsider,
            optionalCommandlineParams->minMapqOfAnchorRead, optionalCommandlineParams->maxMapqOfInrepeatRead);

        return ProgramParameters(paths, readLength, heuristics);
    }

    return boost::none;
}