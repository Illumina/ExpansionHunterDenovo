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

#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ReferenceContigInfo.hh"

class GenomicRegion
{
public:
    friend std::ostream& operator<<(std::ostream& ostrm, const GenomicRegion& region);

    GenomicRegion(int contigId, int64_t start, int64_t end);

    bool operator<(const GenomicRegion& other) const;

    bool overlaps(const GenomicRegion& other) const;
    int64_t distance(const GenomicRegion& other) const;

    int contigId() const { return contigId_; }
    int64_t start() const { return start_; }
    int64_t end() const { return end_; }

    void setContigId(int contigId) { contigId_ = contigId; }
    void setStart(int64_t start) { start_ = start; }
    void setEnd(int64_t end) { end_ = end; }
    std::string asString(const ReferenceContigInfo& contigInfo) const;
    bool operator==(const GenomicRegion& other) const { return equalTo(other); }

protected:
    bool equalTo(const GenomicRegion& other) const
    {
        return contigId_ == other.contigId_ && start_ == other.start_ && end_ == other.end_;
    }

private:
    int contigId_; // value -1 means that the region is unaligned
    int64_t start_;
    int64_t end_;
};

template <typename F> class RegionWithFeature : public GenomicRegion
{
public:
    RegionWithFeature(int contigId, int64_t start, int64_t end, F feature)
        : GenomicRegion(contigId, start, end)
        , feature_(std::move(feature))
    {
    }

    F& feature() { return feature_; }
    const F& feature() const { return feature_; }

    bool operator==(const RegionWithFeature<F>& other) const
    {
        return feature_ == other.feature_ && GenomicRegion::equalTo(other);
    }

private:
    F feature_;
};

class CountFeature
{
public:
    CountFeature(int value)
        : value_(value)
    {
    }

    int value() const { return value_; }
    void combine(CountFeature other) { value_ += other.value_; }

    bool operator==(const CountFeature& other) const { return value_ == other.value_; }

private:
    int value_;
};

class SampleCountFeature
{
public:
    SampleCountFeature(std::unordered_map<std::string, int> value)
        : value_(std::move(value))
    {
    }

    const std::unordered_map<std::string, int>& value() const { return value_; }
    void combine(const SampleCountFeature& other);

    bool operator==(const SampleCountFeature& other) const { return value_ == other.value_; }

private:
    std::unordered_map<std::string, int> value_;
};

using RegionWithCount = RegionWithFeature<CountFeature>;
using RegionWithSampleCount = RegionWithFeature<SampleCountFeature>;

RegionWithCount createCountableRegion(int contigId, int64_t start, int64_t end);

template <typename F> void sortAndMerge(std::vector<RegionWithFeature<F>>& regions, int maxMergeDistance = 500)
{
    using Region = RegionWithFeature<F>;
    std::sort(regions.begin(), regions.end());

    std::unique_ptr<Region> mergedRegionPtr;
    std::vector<Region> mergedRegions;

    for (const auto& region : regions)
    {
        if (!mergedRegionPtr)
        {
            mergedRegionPtr.reset(new Region(region));
        }
        else if (mergedRegionPtr->distance(region) <= maxMergeDistance)
        {
            auto maxEnd = std::max<int64_t>(mergedRegionPtr->end(), region.end());
            mergedRegionPtr->setEnd(maxEnd);
            mergedRegionPtr->feature().combine(region.feature());
        }
        else
        {
            mergedRegions.push_back(std::move(*mergedRegionPtr));
            mergedRegionPtr.reset(new Region(region));
        }
    }

    if (mergedRegionPtr)
    {
        mergedRegions.push_back(std::move(*mergedRegionPtr));
    }

    regions = mergedRegions;
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& region);
std::ostream& operator<<(std::ostream& out, const CountFeature& feature);
std::ostream& operator<<(std::ostream& out, const SampleCountFeature& feature);

template <typename F> std::ostream& operator<<(std::ostream& out, const RegionWithFeature<F>& region)
{
    const GenomicRegion& baseRegion = region;
    out << baseRegion << "(" << region.feature() << ")";
    return out;
}

GenomicRegion decode(const ReferenceContigInfo& contigInfo, const std::string& encoding);
