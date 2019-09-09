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

#include "merge/MultisampleProfile.hh"

using std::vector;

void normalize(MultisampleAnchoredIrrProfile& profile)
{
    for (auto& motifAndSampleCounts : profile)
    {
        sortAndMerge(motifAndSampleCounts.second);
    }
}

void add(
    const SampleId& sampleId, const Motif& motif, const GenomicRegion& region, int numAnchoredIrrs,
    MultisampleAnchoredIrrProfile& profile)
{
    SampleCountFeature sampleCount({ { sampleId, numAnchoredIrrs } });
    RegionWithSampleCount regionWithSampleCount(region.contigId(), region.start(), region.end(), sampleCount);
    profile[motif].push_back(regionWithSampleCount);
}
