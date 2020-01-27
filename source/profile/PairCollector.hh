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

#include <fstream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "reads/Read.hh"
#include "region/GenomicRegion.hh"

enum class ReadType
{
    kIrrRead,
    kAnchorRead,
    kOtherRead
};

enum class PairType
{
    kIrrAnchorPair,
    kIrrIrrPair,
    kOtherPair
};

class ReadCache
{
public:
    bool isReadCached(const Read& read);
    ReadType typeOfRead(const Read& read);
    void eraseRead(const Read& read);
    RegionWithCount extractRegionOfIrrOrAnchor(const Read& read);
    std::string extractUnitOfIrr(const Read& read);
    void cacheAnchorRead(const Read& read);
    void cacheInrepeatRead(const Read& read, const std::string& unit);
    void cacheOtherRead(const Read& read);
    std::string printStats();

private:
    std::unordered_map<std::string, ReadType> readTypes_;
    std::unordered_map<std::string, RegionWithCount> irrAndAnchorLocations_;
    std::unordered_map<std::string, std::string> irrUnits_;
};

class PairCollector
{
public:
    PairCollector(ReferenceContigInfo contigInfo)
        : contigInfo_(std::move(contigInfo))
    {
    }
    ~PairCollector();
    void addAnchor(const Read& read);
    void addIrr(const Read& read, const std::string& unit);
    void addOtherRead(const Read& read);
    std::string PrintStats();
    const std::unordered_map<std::string, std::vector<RegionWithCount>>& anchorRegions() { return anchorRegions_; }
    const std::unordered_map<std::string, std::vector<RegionWithCount>>& irrRegions() { return irrRegions_; };

    void enableReadLogging(const std::string& pathToReadLog);

private:
    void logIrrPair(
        const std::string& fragName, const GenomicRegion& readRegion, const std::string& readUnit,
        const GenomicRegion& mateRegion, const std::string& mateUnit);

    void logAnchoredIrr(
        const std::string& fragName, const std::string& unit, const GenomicRegion& irrRegion,
        const GenomicRegion& anchorRegion);

    ReferenceContigInfo contigInfo_;
    ReadCache unparedCache_;
    // Regions containing anchors and IRRs.
    std::unordered_map<std::string, std::vector<RegionWithCount>> anchorRegions_;
    std::unordered_map<std::string, std::vector<RegionWithCount>> irrRegions_;

    std::unique_ptr<std::ofstream> logStream_;
};
