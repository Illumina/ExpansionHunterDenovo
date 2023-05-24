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

#include "SampleRunStats.hh"

#include <algorithm>
#include <sstream>
#include <stdexcept>

using boost::optional;
using std::vector;

bool SampleRunStats::operator==(const SampleRunStats& other) const
{
    return meanReadLength_ == other.meanReadLength_ && depth_ == other.depth_;
}

std::ostream& operator<<(std::ostream& out, const SampleRunStats& stats)
{
    out << "SampleRunStats(meanReadLength=" << stats.meanReadLength() << ", depth=" << stats.depth() << ")";
    return out;
}

SampleRunStatsCalculator::SampleRunStatsCalculator(ReferenceContigInfo contigInfo)
    : contigInfo_(std::move(contigInfo))
    , totalReadCount(0)
    , sumOfReadLengths(0)
{
}

void SampleRunStatsCalculator::inspect(int contigId, int readLength)
{
    ++totalReadCount;
    sumOfReadLengths += readLength;
    static int kNumContigsToConsider = 22;
    if (contigId <= kNumContigsToConsider)
    {
        ++contigIdToReadCount[contigId];
    }
}

static double median(vector<double> numbers)
{
    if (numbers.empty())
    {
        throw std::logic_error("Median of an empty array is undefined");
    }

    std::sort(numbers.begin(), numbers.end());
    if (numbers.size() % 2 == 0)
    {
        double first = numbers[numbers.size() / 2 - 1];
        double second = numbers[numbers.size() / 2];
        return ((first + second) / 2);
    }
    else
    {
        return numbers[numbers.size() / 2];
    }
}

optional<SampleRunStats> SampleRunStatsCalculator::estimate() const
{
    if (totalReadCount == 0)
    {
        return boost::none;
    }

    const auto meanReadLength = static_cast<int>(sumOfReadLengths / totalReadCount);

    vector<double> meanDepths;
    for (const auto& contigIdAndReadCount : contigIdToReadCount)
    {
        const auto contigId = contigIdAndReadCount.first;
        // Skip unaligned reads
        if (contigId == -1)
        {
            continue;
        }

        const auto readCount = contigIdAndReadCount.second;
        const auto contigLength = contigInfo_.getContigSize(contigId);
        const auto medianDepth = readCount * meanReadLength / static_cast<double>(contigLength);

        meanDepths.push_back(medianDepth);
    }

    const double medianDepth = median(meanDepths);

    return SampleRunStats(meanReadLength, medianDepth);
}
