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

#include "region/GenomicRegion.hh"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using std::ostream;
using std::string;
using std::vector;

GenomicRegion::GenomicRegion(int contigId, int64_t start, int64_t end)
    : contigId_(contigId)
    , start_(start)
    , end_(end)
{
    if (contigId_ < -1)
    {
        std::ostringstream out;
        out << *this;
        throw std::logic_error("Region " + out.str() + " is invalid");
    }
}

bool GenomicRegion::operator<(const GenomicRegion& other) const
{
    if (contigId_ != other.contigId_)
    {
        return contigId_ < other.contigId_;
    }

    if (start_ != other.start_)
    {
        return start_ < other.start_;
    }

    return end_ < other.end_;
}

bool GenomicRegion::overlaps(const GenomicRegion& other) const
{
    if (contigId_ != other.contigId_)
    {
        return false;
    }

    // All unaligned regions assumed to overlap.
    if (contigId_ == -1)
    {
        return true;
    }

    const int64_t left_bound = start_ > other.start_ ? start_ : other.start_;
    const int64_t right_bound = end_ < other.end_ ? end_ : other.end_;

    return left_bound <= right_bound;
}

string GenomicRegion::asString(const ReferenceContigInfo& contigInfo) const
{
    if (contigId_ == -1)
    {
        return "unaligned";
    }

    const string& contigName = contigInfo.getContigName(contigId_);
    std::ostringstream out;
    out << contigName << ':' << start_ << '-' << end_;
    return out.str();
}

int64_t GenomicRegion::distance(const GenomicRegion& other) const
{
    if (contigId_ != other.contigId_)
    {
        return std::numeric_limits<int64_t>::max();
    }

    // All unaligned regions assumed to overlap.
    if (contigId_ == -1)
    {
        return 0;
    }

    if (end_ < other.start_)
    {
        return other.start_ - end_;
    }

    if (other.end_ < start_)
    {
        return start_ - other.end_;
    }

    return 0;
}

void SampleCountFeature::combine(const SampleCountFeature& other)
{
    for (const auto& sampleCount : other.value())
    {
        value_[sampleCount.first] += sampleCount.second;
    }
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& region)
{
    out << region.contigId_ << ":" << region.start_ << "-" << region.end_;
    return out;
}

std::ostream& operator<<(std::ostream& out, const CountFeature& feature)
{
    out << feature.value();
    return out;
}

std::ostream& operator<<(std::ostream& out, const SampleCountFeature& feature)
{
    bool firstElement = true;
    for (const auto& sampleCount : feature.value())
    {
        if (firstElement)
        {
            firstElement = false;
        }
        else
        {
            out << ", ";
        }

        out << "{" << sampleCount.first << ", " << sampleCount.second << "}";
    }
    return out;
}

RegionWithCount createCountableRegion(int contigId, int64_t start, int64_t end)
{
    return { contigId, start, end, CountFeature(1) };
}

GenomicRegion decode(const ReferenceContigInfo& contigInfo, const string& encoding)
{
    if (encoding == "unaligned")
    {
        return { -1, 0, 0 };
    }

    auto colonIndex = encoding.find_last_of(':');
    if (colonIndex == string::npos || colonIndex == 0 || colonIndex + 1 == encoding.size())
    {
        throw std::logic_error("Unexpected range format: " + encoding);
    }

    string contig = encoding.substr(0, colonIndex);
    string intervalEncoding = encoding.substr(colonIndex + 1);

    auto numDashes = std::count(intervalEncoding.begin(), intervalEncoding.end(), '-');
    if (numDashes != 1)
    {
        throw std::logic_error("Unexpected range format: " + encoding);
    }

    auto dashIndex = intervalEncoding.find('-');
    if (dashIndex == 0 || dashIndex + 1 == intervalEncoding.size())
    {
        throw std::logic_error("Unexpected range format: " + encoding);
    }

    int contigIndex = contigInfo.getContigId(contig);
    int64_t start = std::stoi(intervalEncoding.substr(0, dashIndex));
    int64_t end = std::stoi(intervalEncoding.substr(dashIndex + 1));

    return { contigIndex, start, end };
}
