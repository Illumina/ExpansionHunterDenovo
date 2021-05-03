//
// ExpansionHunter Denovo
// Copyright 2016-2020 Illumina, Inc.
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

#include "outlier/OutlierWorkflow.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "thirdparty/nlohmann_json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "common/Manifest.hh"
#include "io/Reference.hh"
#include "outlier/ZScoreAnalysis.hh"
#include "region/GenomicRegion.hh"

using Json = nlohmann::json;
using std::string;
using std::unordered_map;
using std::vector;

struct SampleStats
{
    SampleStats(double depth, int readLength)
        : depth(depth)
        , readLength(readLength)
    {
    }
    double depth;
    int readLength;
};

using StatsBySample = unordered_map<string, SampleStats>;

vector<AnchoredIrrCounts> extractAnchoredIrrCounts(const Reference& reference, const Json& multisampleProfile)
{
    if (!multisampleProfile.contains("Counts"))
    {
        throw std::runtime_error("Malformed multisample profile: Counts are missing");
    }

    vector<AnchoredIrrCounts> irrCountsByRegionAndMotif;
    const auto& counts = multisampleProfile["Counts"];

    for (const auto& motifAndFindings : counts.items())
    {
        const string& motif = motifAndFindings.key();
        const auto& findings = motifAndFindings.value();
        if (!findings.contains("RegionsWithIrrAnchors"))
        {
            continue;
        }

        for (const auto& regionAndCountsBySample : findings["RegionsWithIrrAnchors"].items())
        {
            const string& regionEncoding = regionAndCountsBySample.key();
            const auto& countsBySample = regionAndCountsBySample.value();

            GenomicRegion region = decode(reference.contigInfo(), regionEncoding);
            AnchoredIrrCounts anchoredIrrCounts(region, motif);
            for (const auto& sampleAndCount : countsBySample.items())
            {
                const string& sample = sampleAndCount.key();
                double count = sampleAndCount.value();
                anchoredIrrCounts.countBySample.emplace(sample, count);
            }
            irrCountsByRegionAndMotif.emplace_back(anchoredIrrCounts);
        }
    }

    return irrCountsByRegionAndMotif;
}

StatsBySample extractSampleStats(const Json& multisampleProfile)
{
    unordered_map<string, SampleStats> statsBySample;
    if (!multisampleProfile.contains("Parameters"))
    {
        throw std::runtime_error("Malformed multisample profile: Sample parameters are missing");
    }

    const auto& params = multisampleProfile["Parameters"];

    if (!params.contains("Depths"))
    {
        throw std::runtime_error("Malformed multisample profile: Sample parameters are missing depths");
    }
    const auto& depths = params["Depths"];
    for (const auto& sampleAndDepth : depths.items())
    {
        const string& sample = sampleAndDepth.key();
        double depth = sampleAndDepth.value();
        statsBySample.emplace(sample, SampleStats(depth, -1));
    }

    if (!params.contains("ReadLengths"))
    {
        throw std::runtime_error("Malformed multisample profile: Sample parameters are missing read lengths");
    }

    const auto& readLengths = params["ReadLengths"];
    for (const auto& sampleAndReadLength : readLengths.items())
    {
        const string& sample = sampleAndReadLength.key();
        int readLength = sampleAndReadLength.value();
        auto record = statsBySample.find(sample);
        if (record == statsBySample.end())
        {
            throw std::runtime_error("Malformed multisample profile: Missing read length for " + sample);
        }
        record->second.readLength = readLength;
    }

    return statsBySample;
}

void depthNormalize(
    const StatsBySample& statsBySample, double targetDepth, vector<AnchoredIrrCounts>& anchoredIrrCountTable)
{
    for (auto& anchoredIrrCount : anchoredIrrCountTable)
    {
        for (auto& sampleAndCount : anchoredIrrCount.countBySample)
        {
            const auto& sample = sampleAndCount.first;
            double& count = sampleAndCount.second;
            auto record = statsBySample.find(sample);
            if (record == statsBySample.end())
            {
                throw std::runtime_error("Missing stats for sample " + sample);
            }

            const double actualDepth = record->second.depth;
            count = targetDepth * count / actualDepth;
        }
    }
}

void performLocusAnalysis(
    const Reference& genome, const Manifest& manifest, const vector<AnchoredIrrCounts>& irrCountsByRegionAndMotif,
    const string& outputPath)
{
    std::ofstream outputFile(outputPath);
    if (!outputFile)
    {
        throw std::logic_error("Unable to write to " + outputPath);
    }
    outputFile << std::fixed << std::setprecision(2);
    outputFile << "contig\tstart\tend\tmotif\ttop_case_zscore\thigh_case_counts\tnon_zero_counts" << std::endl;
    for (const auto& irrCounts : irrCountsByRegionAndMotif)
    {
        const auto& region = irrCounts.region;
        const auto& contig = genome.contigInfo().getContigName(region.contigId());

        auto results = analyzeZScores(manifest, irrCounts);

        if (results.casesWithHighCounts.empty())
        {
            continue;
        }

        outputFile << contig << "\t" << region.start() << "\t" << region.end() << "\t" << irrCounts.motif << "\t";
        outputFile << results.topZScore << "\t";
        for (int sampleIndex = 0; sampleIndex != results.casesWithHighCounts.size(); ++sampleIndex)
        {
            const auto& sample = results.casesWithHighCounts[sampleIndex];
            outputFile << sample << ":" << irrCounts.countBySample.find(sample)->second;
            if (sampleIndex + 1 != results.casesWithHighCounts.size())
            {
                outputFile << ",";
            }
        }
        outputFile << "\t";

        for (auto record = irrCounts.countBySample.begin(); record != irrCounts.countBySample.end(); ++record)
        {
            outputFile << record->second;
            if (std::next(record) != irrCounts.countBySample.end())
            {
                outputFile << ",";
            }
        }

        outputFile << std::endl;
    }
}

int runOutlierWorkflow(const OutlierWorkflowParameters& parameters)
{
    Reference reference(parameters.pathToReference());

    std::ifstream multisampleProfileFile(parameters.pathToMultisampleProfile());
    if (!multisampleProfileFile)
    {
        throw std::runtime_error("Unable to read " + parameters.pathToMultisampleProfile());
    }

    Json multisampleProfile;
    multisampleProfileFile >> multisampleProfile;

    auto anchoredIrrCountsByRegion = extractAnchoredIrrCounts(reference, multisampleProfile);
    spdlog::info("Loaded {} regions", anchoredIrrCountsByRegion.size());

    spdlog::info("Normalizing counts");
    auto sampleStats = extractSampleStats(multisampleProfile);

    for (const auto& regionCounts : anchoredIrrCountsByRegion)
    {
        std::cerr << regionCounts.region << " " << regionCounts.motif << std::endl;
        for (const auto& sampleAndCount : regionCounts.countBySample)
        {
            std::cerr << " " << sampleAndCount.first << " " << sampleAndCount.second << std::endl;
        }
    }

    auto statsBySample = extractSampleStats(multisampleProfile);

    for (const auto& sampleAndStats : statsBySample)
    {
        std::cerr << sampleAndStats.first << "\t" << sampleAndStats.second.depth << "\t"
                  << sampleAndStats.second.readLength << std::endl;
    }

    const double targetDepth = 40;
    depthNormalize(statsBySample, targetDepth, anchoredIrrCountsByRegion);

    for (const auto& regionCounts : anchoredIrrCountsByRegion)
    {
        std::cerr << regionCounts.region << " " << regionCounts.motif << std::endl;
        for (const auto& sampleAndCount : regionCounts.countBySample)
        {
            std::cerr << " " << sampleAndCount.first << " " << sampleAndCount.second << std::endl;
        }
    }

    vector<string> orderedSamples;
    auto manifest = loadManifest(parameters.pathToManifest(), orderedSamples);
    performLocusAnalysis(reference, manifest, anchoredIrrCountsByRegion, parameters.pathToLocusResults());

    return 0;
}
