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
using SampleIdToSampleParameters = unordered_map<string, SampleParameters>;

enum class SampleStatus
{
    kCase,
    kControl
};

SampleStatus decodeSampleStatus(const string& encoding)
{
    if (encoding == "case")
    {
        return SampleStatus::kCase;
    }
    else if (encoding == "control")
    {
        return SampleStatus::kControl;
    }
    else
    {
        throw std::runtime_error(encoding + " is not a valid sample status");
    }
}

struct ManifestEntry
{
    ManifestEntry(string sample, SampleStatus status, string path)
        : sample(std::move(sample))
        , status(status)
        , path(std::move(path))
    {
    }

    ManifestEntry(string sample, const string& statusEncoding, string path)
        : sample(std::move(sample))
        , status(decodeSampleStatus(statusEncoding))
        , path(std::move(path))
    {
    }

    string sample;
    SampleStatus status;
    string path;
};

using Manifest = vector<ManifestEntry>;

Manifest loadManifest(const string& path)
{
    Manifest manifest;

    std::ifstream manifestFile(path);
    if (!manifestFile)
    {
        throw std::runtime_error("Unable to load manifest from " + path);
    }

    string line;
    while (std::getline(manifestFile, line))
    {
        std::istringstream decoder(line);
        string sample;
        string statusEncoding;
        string path;
        if (!(decoder >> sample >> statusEncoding >> path))
        {
            throw std::runtime_error("Unable to decode manifest line " + line);
        }
        manifest.emplace_back(sample, statusEncoding, path);
    }

    return manifest;
}

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
    const ManifestEntry& sampleInfo, const ReferenceContigInfo& contigInfo,
    MultisampleAnchoredIrrProfile& anchoredIrrProfile, MultisampleIrrPairProfile& pairedIrrProfile,
    SampleIdToSampleParameters& parametersForSamples)
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
        else
        {
            const string& motif = record.key();
            loadAnchorInfo(contigInfo, sampleInfo.sample, motif, record.value(), anchoredIrrProfile);
            loadPairedIrrProfile(sampleInfo.sample, motif, record.value(), pairedIrrProfile);
        }
    }

    if (readLength == 0)
    {
        throw std::runtime_error("Read length appears to be unset for " + sampleInfo.sample);
    }

    if (depth == -1)
    {
        throw std::runtime_error("Depth appears to be unset for " + sampleInfo.sample);
    }

    parametersForSamples.emplace(sampleInfo.sample, SampleParameters(readLength, depth));
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

    Manifest manifest = loadManifest(parameters.pathToManifest());
    spdlog::info("Loaded manifest describing {} samples", manifest.size());

    MultisampleAnchoredIrrProfile anchoredIrrProfile;
    MultisampleIrrPairProfile irrPairProfile;
    SampleIdToSampleParameters parametersForSamples;

    int sampleCount = 0;
    int kNormalizationStride = 50;
    for (const auto& sampleInfo : manifest)
    {
        spdlog::info("Loading STR profile of {}", sampleInfo.sample);
        loadSampleProfile(sampleInfo, contigInfo, anchoredIrrProfile, irrPairProfile, parametersForSamples);
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