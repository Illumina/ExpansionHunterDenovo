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

#include "ReferenceContigInfo.hh"

#include <memory>
#include <stdexcept>

using std::pair;
using std::string;
using std::to_string;
using std::vector;

// Removes "chr" prefix from contig names that contain it; adds it to contigs that don't
static string generateAlternativeContigName(const string& originalName)
{
    if (originalName.length() > 3 && originalName.substr(0, 3) == "chr")
    {
        return originalName.substr(3);
    }
    else
    {
        return "chr" + originalName;
    }
}

ReferenceContigInfo::ReferenceContigInfo(vector<pair<string, int64_t>> namesAndSizes)
    : namesAndSizes_(std::move(namesAndSizes))
{
    for (int index = 0; index != static_cast<int>(namesAndSizes_.size()); ++index)
    {
        const auto& contigName = namesAndSizes_[index].first;
        nameToIndex_.emplace(std::make_pair(contigName, index));
    }
}

const std::string& ReferenceContigInfo::getContigName(int contigIndex) const
{
    assertValidIndex(contigIndex);
    return namesAndSizes_[contigIndex].first;
}

int64_t ReferenceContigInfo::getContigSize(int contigIndex) const
{
    assertValidIndex(contigIndex);
    return namesAndSizes_[contigIndex].second;
}

int ReferenceContigInfo::getContigId(const std::string& contigName) const
{
    auto entry = nameToIndex_.find(contigName);
    if (entry == nameToIndex_.end())
    {
        entry = nameToIndex_.find(generateAlternativeContigName(contigName));
    }

    if (entry == nameToIndex_.end())
    {
        throw std::logic_error("Invalid contig name " + contigName);
    }

    return entry->second;
}

void ReferenceContigInfo::assertValidIndex(int contigIndex) const
{
    if (contigIndex < 0 || contigIndex >= static_cast<int>(namesAndSizes_.size()))
    {
        throw std::logic_error("Invalid contig index " + to_string(contigIndex));
    }
}

std::ostream& operator<<(std::ostream& out, const ReferenceContigInfo& contigInfo)
{
    for (int contigIndex = 0; contigIndex != contigInfo.numContigs(); ++contigIndex)
    {
        const auto& contigName = contigInfo.getContigName(contigIndex);
        out << contigName << " -> " << contigIndex << std::endl;
    }

    return out;
}
