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

#include "merge/MergeWorkflow.hh"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "thirdparty/nlohmann_json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "MergeParameters.hh"
#include "common/Manifest.hh"
#include "io/Reference.hh"
#include "merge/MultisampleProfile.hh"

using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using Json = nlohmann::json;

struct SampleParameters
{
    SampleParameters(int readLength, double depth)
        : readLength(readLength)
        , depth(depth)
    {
    }
    int readLength;
    double depth;
};
using SampleIdToSampleParameters = std::unordered_map<std::string, SampleParameters>;

void loadAnchorInfo(
    const ReferenceContigInfo& contigInfo, const string& sampleId, const string& motif, const Json& record,
    MultisampleAnchoredIrrProfile& anchoredIrrProfile)
{
    if (record.find("RegionsWithIrrAnchors") == record.end())
    {
        return;
    }

    for (const auto& regionAndCount : record["RegionsWithIrrAnchors"].items())
    {
        GenomicRegion region = decode(contigInfo, regionAndCount.key());
        SampleCountFeature sampleCount({ { sampleId, regionAndCount.value() } });
        RegionWithSampleCount regionWithSampleCount(region.contigId(), region.start(), region.end(), sampleCount);
        anchoredIrrProfile[motif].push_back(regionWithSampleCount);
    }
}

void loadPairedIrrProfile(
    const string& sampleId, const string& motif, const Json& record, MultisampleIrrPairProfile& pairedIrrProfile)
{
    if (record.find("IrrPairCount") == record.end() || record["IrrPairCount"] == 0)
    {
        return;
    }

    pairedIrrProfile[motif].emplace(sampleId, record["IrrPairCount"]);
}

void loadSampleProfile(
    const std::string& sample, const ManifestEntry& sampleInfo, const ReferenceContigInfo& contigInfo,
    MultisampleAnchoredIrrProfile& anchoredIrrProfile, MultisampleIrrPairProfile& pairedIrrProfile,
    SampleIdToSampleParameters& parametersForSamples, int shortestUnit, int longestUnit)
{
    std::ifstream profileFile(sampleInfo.path);
    if (!profileFile)
    {
        throw std::runtime_error("Unable to read " + sampleInfo.path);
    }

    Json profileJson;
    profileFile >> profileJson;

    int readLength = 0;
    double depth = -1;
    for (const auto& record : profileJson.items())
    {
        if (record.key() == "ReadLength")
        {
            readLength = record.value();
        }
        else if (record.key() == "Depth")
        {
            depth = record.value();
        }
        else if (shortestUnit <= record.key().length() && record.key().length() <= longestUnit)
        {
            const string& motif = record.key();
            loadAnchorInfo(contigInfo, sample, motif, record.value(), anchoredIrrProfile);
            loadPairedIrrProfile(sample, motif, record.value(), pairedIrrProfile);
        }
    }

    if (readLength == 0)
    {
        throw std::runtime_error("Read length appears to be unset for " + sample);
    }

    if (depth == -1)
    {
        throw std::runtime_error("Depth appears to be unset for " + sample);
    }

    parametersForSamples.emplace(sample, SampleParameters(readLength, depth));
}

void writeMultisampleProfile(
    const ReferenceContigInfo& contigInfo, const string& outputPath,
    const MultisampleAnchoredIrrProfile& anchoredIrrProfile, const MultisampleIrrPairProfile& pairedIrrProfile,
    const SampleIdToSampleParameters& parametersForSamples)
{
    Json countsRecord;
    for (const auto& motifAndRecord : pairedIrrProfile)
    {
        const string& motif = motifAndRecord.first;
        const auto& record = motifAndRecord.second;
        for (const auto& sampleIdAndCount : record)
        {
            const string& sampleId = sampleIdAndCount.first;
            int count = sampleIdAndCount.second;
            countsRecord[motif]["IrrPairCounts"][sampleId] = count;
        }
    }

    for (const auto& motifAndRecord : anchoredIrrProfile)
    {
        const string& motif = motifAndRecord.first;
        const auto& record = motifAndRecord.second;
        for (const auto& regionWithSampleCounts : record)
        {
            const auto& regionEncoding = regionWithSampleCounts.asString(contigInfo);
            for (const auto& sampleIdAndCount : regionWithSampleCounts.feature().value())
            {
                const auto& sampleId = sampleIdAndCount.first;
                int count = sampleIdAndCount.second;
                countsRecord[motif]["RegionsWithIrrAnchors"][regionEncoding][sampleId] = count;
            }
        }
    }

    Json parametersRecord;
    for (const auto& sampleIdAndParameters : parametersForSamples)
    {
        const auto& sampleId = sampleIdAndParameters.first;
        const auto& parameters = sampleIdAndParameters.second;
        parametersRecord["ReadLengths"][sampleId] = parameters.readLength;
        parametersRecord["Depths"][sampleId] = parameters.depth;
    }

    Json multisampleProfile;
    multisampleProfile["Counts"] = countsRecord;
    multisampleProfile["Parameters"] = parametersRecord;

    std::ofstream outputFile(outputPath);
    if (!outputFile)
    {
        throw std::logic_error("Unable to write to " + outputPath);
    }
    outputFile << multisampleProfile.dump(4);
}

int runMergeWorkflow(const MergeWorkflowParameters& parameters)
{
    assertValidity(parameters);
    Reference reference(parameters.pathToReference());
    const ReferenceContigInfo& contigInfo = reference.contigInfo();

    vector<string> orderedSamples;
    Manifest manifest = loadManifest(parameters.pathToManifest(), orderedSamples);
    spdlog::info("Loaded manifest describing {} samples", manifest.size());

    MultisampleAnchoredIrrProfile anchoredIrrProfile;
    MultisampleIrrPairProfile irrPairProfile;
    SampleIdToSampleParameters parametersForSamples;

    int sampleCount = 0;
    int kNormalizationStride = 50;
    for (const auto& sample : orderedSamples)
    {
        const auto& sampleInfo = manifest.find(sample)->second;
        spdlog::info("Loading STR profile of {}", sample);
        loadSampleProfile(
            sample, sampleInfo, contigInfo, anchoredIrrProfile, irrPairProfile, parametersForSamples,
            parameters.shortestUnitToConsider(), parameters.longestUnitToConsider());
        sampleCount++;

        if (sampleCount % kNormalizationStride == 0)
        {
            spdlog::info("Normalizing after loading sample #{}", sampleCount);
            normalize(anchoredIrrProfile);
        }
    }
    normalize(anchoredIrrProfile);

    writeMultisampleProfile(
        contigInfo, parameters.pathToMultisampleProfile(), anchoredIrrProfile, irrPairProfile, parametersForSamples);

    spdlog::info("Done");
    return 0;
}