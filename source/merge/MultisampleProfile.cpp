//
// ExpansionHunter Denovo
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Michael Eberle <meberle@illumina.com>
//
// Licensed under the PolyForm Strict License 1.0.0
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      https://polyformproject.org/licenses/strict/1.0.0
//
// As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.
// See the License for the specific language governing permissions and
// limitations under the License.

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
