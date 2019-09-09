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

#include "profile/ProfileWorkflow.hh"

#include <fstream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "thirdparty/nlohmann_json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "io/HtsFileStreamer.hh"
#include "profile/PairCollector.hh"
#include "profile/ReadClassification.hh"
#include "profile/SampleRunStats.hh"

using std::set;
using std::string;
using std::unordered_map;
using std::vector;

int runProfileWorkflow(const ProfileWorkflowParameters& parameters)
{
    assertValidity(parameters);
    spdlog::info("File with reads: {}", parameters.pathToReads());

    PairCollector pairCollector;
    HtsFileStreamer readStreamer(parameters.pathToReads(), parameters.pathToReference());

    const ReferenceContigInfo& referenceContigInfo = readStreamer.contigInfo();
    SampleRunStatsCalculator statsCalculator(referenceContigInfo);

    while (readStreamer.trySeekingToNextPrimaryAlignment())
    {
        statsCalculator.inspect(readStreamer.currentReadContigId(), readStreamer.currentReadLength());

        Read read = readStreamer.decodeRead();

        string repeatUnit;
        const ReadType readType
            = classifyRead(parameters.maxMapqOfInrepeatRead(), parameters.minMapqOfAnchorRead(), read, repeatUnit);
        if (readType == ReadType::kIrrRead)
        {
            pairCollector.addIrr(read, repeatUnit);
        }
        else if (readType == ReadType::kAnchorRead)
        {
            pairCollector.addAnchor(read);
        }
        else
        {
            pairCollector.addOtherRead(read);
        }
    }

    const auto stats = statsCalculator.estimate();
    nlohmann::json output;
    output["ReadLength"] = stats->meanReadLength();
    output["Depth"] = stats->depth();

    const auto irrAnchorRegions = pairCollector.anchorRegions();
    const auto irrRegions = pairCollector.irrRegions();

    set<string> units;
    for (const auto& kv : irrRegions)
    {
        const string unit = kv.first;
        if (parameters.shortestUnitToConsider() <= unit.length() && unit.length() <= parameters.longestUnitToConsider())
        {
            units.insert(unit);
        }
    }

    for (const auto& unit : units)
    {

        output[unit]["RepeatUnit"] = unit;

        size_t ancIrrCount = 0;
        bool foundAncIrrs = irrAnchorRegions.find(unit) != irrAnchorRegions.end();
        if (foundAncIrrs)
        {
            ancIrrCount = irrAnchorRegions.at(unit).size();
        }

        size_t irrPairCount = (irrRegions.at(unit).size() - ancIrrCount) / 2;

        output[unit]["AnchoredIrrCount"] = ancIrrCount;
        output[unit]["IrrPairCount"] = irrPairCount;

        if (foundAncIrrs)
        {
            vector<RegionWithCount> mergedRegionsWithAnchors = irrAnchorRegions.at(unit);
            sortAndMerge(mergedRegionsWithAnchors);

            for (const auto& region : mergedRegionsWithAnchors)
            {
                const string regionEncoding = region.asString(referenceContigInfo);
                output[unit]["RegionsWithIrrAnchors"][regionEncoding] = region.feature().value();
            }
        }
    }

    const string& jsonPath = parameters.profilePath();
    std::ofstream jsonStream;
    jsonStream.open(jsonPath.c_str());

    if (!jsonStream.is_open())
    {
        throw std::runtime_error(
            "Failed to open output JSON file " + jsonPath + " for writing (" + strerror(errno) + ")");
    }

    jsonStream << output.dump(4) << std::endl;

    return 0;
}