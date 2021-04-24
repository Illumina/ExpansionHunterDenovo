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
#include <iostream>
#include <stdexcept>

#include "thirdparty/nlohmann_json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "io/Reference.hh"
#include "region/GenomicRegion.hh"

using Json = nlohmann::json;
using std::string;
using std::unordered_map;
using std::vector;

struct AnchoredIrrCounts
{
    AnchoredIrrCounts(GenomicRegion region, string motif)
        : region(region)
        , motif(std::move(motif))
    {
    }
    GenomicRegion region;
    string motif;
    unordered_map<string, int> countBySample;
};

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
                int count = sampleAndCount.value();
                anchoredIrrCounts.countBySample.emplace(sample, count);
            }
            irrCountsByRegionAndMotif.emplace_back(anchoredIrrCounts);
        }
    }

    return irrCountsByRegionAndMotif;
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

    for (const auto& regionCounts : anchoredIrrCountsByRegion)
    {
        std::cerr << regionCounts.region << " " << regionCounts.motif << std::endl;
        for (const auto& sampleAndCount : regionCounts.countBySample)
        {
            std::cerr << " " << sampleAndCount.first << " " << sampleAndCount.second << std::endl;
        }
    }

    // std::cerr << multisampleProfile;

    return 0;
}
