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
#include <iomanip>
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

using RegionsByUnit = std::unordered_map<std::string, std::vector<RegionWithCount>>;

set<string> getTargetRepeatUnits(const RegionsByUnit& irrRegions, Interval targetSizeRange)
{
    set<string> units;
    for (const auto& unitAndRegions : irrRegions)
    {
        const string& unit = unitAndRegions.first;
        if (targetSizeRange.contains(unit.length()))
        {
            units.insert(unit);
        }
    }

    return units;
}

void outputProfile(
    const string& profilePath, const SampleRunStats& sampleStats, const RegionsByUnit& irrAnchorRegions,
    const RegionsByUnit& irrRegions, const set<string>& targetUnits, const ReferenceContigInfo& contigInfo)
{
    nlohmann::json output;
    output["ReadLength"] = sampleStats.meanReadLength();
    output["Depth"] = sampleStats.depth();

    for (const auto& unit : targetUnits)
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
                const string regionEncoding = region.asString(contigInfo);
                output[unit]["RegionsWithIrrAnchors"][regionEncoding] = region.feature().value();
            }
        }
    }

    std::ofstream profileStream;
    profileStream.open(profilePath.c_str());

    if (!profileStream.is_open())
    {
        throw std::runtime_error(
            "Failed to open output JSON file " + profilePath + " for writing (" + strerror(errno) + ")");
    }

    profileStream << output.dump(4) << std::endl;
    profileStream.close();
}

int getNumUnitsSpanned(double sampleDepth, int readLength, int unitLen, int numIrrs)
{
    const double alleleDepth = sampleDepth / 2;
    const double numStarts = (numIrrs * readLength) / alleleDepth;
    const double lengthInBp = readLength + numStarts;
    const double lengthInUnits = lengthInBp / unitLen;
    return static_cast<int>(lengthInUnits);
}

double depthNormalize(double sampleDepth, int numIrrs, double targetDepth = 30.0)
{
    return (numIrrs * targetDepth) / sampleDepth;
}

void outputLocusTable(
    const string& tablePath, const SampleRunStats& sampleStats, const RegionsByUnit& irrAnchorRegions,
    const set<string>& targetUnits, const ReferenceContigInfo& contigInfo)
{
    std::ofstream tableStream;
    tableStream.open(tablePath.c_str());

    if (!tableStream.is_open())
    {
        throw std::runtime_error("Failed to open output file " + tablePath + " for writing (" + strerror(errno) + ")");
    }

    tableStream << std::setprecision(2) << std::fixed;
    tableStream << "contig\tstart\tend\tmotif\tnum_anc_irrs\tnorm_num_anc_irrs\thet_str_size" << std::endl;
    for (const auto& unit : targetUnits)
    {
        if (irrAnchorRegions.find(unit) == irrAnchorRegions.end())
        {
            continue;
        }

        vector<RegionWithCount> mergedRegionsWithAnchors = irrAnchorRegions.at(unit);
        sortAndMerge(mergedRegionsWithAnchors);

        for (const auto& region : mergedRegionsWithAnchors)
        {
            if (region.contigId() == -1)
            {
                continue;
            }

            const string& contigName = contigInfo.getContigName(region.contigId());
            const int numIrrs = region.feature().value();
            const double normNumIrrs = depthNormalize(sampleStats.depth(), numIrrs);

            const int numUnitsSpanned
                = getNumUnitsSpanned(sampleStats.depth(), sampleStats.meanReadLength(), unit.length(), numIrrs);
            tableStream << contigName << "\t" << region.start() << "\t" << region.end() << "\t" << unit << "\t"
                        << numIrrs << "\t" << normNumIrrs << "\t" << numUnitsSpanned << std::endl;
        }
    }

    tableStream.close();
}

void outputMotifTable(
    const string& tablePath, const SampleRunStats& sampleStats, const RegionsByUnit& irrAnchorRegions,
    const RegionsByUnit& irrRegions, const set<string>& targetUnits, const ReferenceContigInfo& contigInfo)
{
    std::ofstream tableStream;
    tableStream.open(tablePath.c_str());

    if (!tableStream.is_open())
    {
        throw std::runtime_error("Failed to open output file " + tablePath + " for writing (" + strerror(errno) + ")");
    }

    tableStream << std::setprecision(2) << std::fixed;
    tableStream << "motif\tnum_paired_irrs\tnorm_num_paired_irrs" << std::endl;
    for (const auto& unit : targetUnits)
    {
        int ancIrrCount = 0;
        if (irrAnchorRegions.find(unit) != irrAnchorRegions.end())
        {
            ancIrrCount = irrAnchorRegions.at(unit).size();
        }

        const int irrPairCount = (static_cast<int>(irrRegions.at(unit).size()) - ancIrrCount) / 2;
        if (irrPairCount == 0)
        {
            continue;
        }

        const double normIrrPairCount = depthNormalize(sampleStats.depth(), irrPairCount);
        tableStream << unit << "\t" << irrPairCount << "\t" << normIrrPairCount << std::endl;
    }
    tableStream.close();
}

int runProfileWorkflow(const ProfileWorkflowParameters& parameters)
{
    assertValidity(parameters);
    spdlog::info("File with reads: {}", parameters.pathToReads());

    HtsFileStreamer readStreamer(parameters.pathToReads(), parameters.pathToReference());

    const ReferenceContigInfo& referenceContigInfo = readStreamer.contigInfo();
    SampleRunStatsCalculator statsCalculator(referenceContigInfo);

    PairCollector pairCollector(referenceContigInfo);
    if (parameters.pathToReadLog())
    {
        pairCollector.enableReadLogging(*parameters.pathToReadLog());
    }

    while (readStreamer.trySeekingToNextPrimaryAlignment())
    {
        statsCalculator.inspect(readStreamer.currentReadContigId(), readStreamer.currentReadLength());

        Read read = readStreamer.decodeRead();

        string motif;
        const ReadType readType = classifyRead(
            parameters.motifSizeRange(), parameters.maxMapqOfInrepeatRead(), parameters.minMapqOfAnchorRead(), read,
            motif);
        if (readType == ReadType::kIrrRead)
        {
            pairCollector.addIrr(read, motif);
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
    assert(stats);

    auto targetUnits = getTargetRepeatUnits(pairCollector.irrRegions(), parameters.motifSizeRange());
    outputProfile(
        parameters.profilePath(), *stats, pairCollector.anchorRegions(), pairCollector.irrRegions(), targetUnits,
        referenceContigInfo);
    outputLocusTable(
        parameters.pathToLocusTable(), *stats, pairCollector.anchorRegions(), targetUnits, referenceContigInfo);
    outputMotifTable(
        parameters.pathToMotifTable(), *stats, pairCollector.anchorRegions(), pairCollector.irrRegions(), targetUnits,
        referenceContigInfo);
    return 0;
}
