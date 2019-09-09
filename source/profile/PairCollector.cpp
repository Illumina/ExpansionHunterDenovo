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

#include "PairCollector.hh"

#include <cassert>

using std::string;
using std::to_string;

bool ReadCache::isReadCached(const Read& read)
{
    const auto it = readTypes_.find(read.name);
    return it != readTypes_.end();
}

ReadType ReadCache::typeOfRead(const Read& read)
{
    const auto it = readTypes_.find(read.name);
    if (it == readTypes_.end())
    {
        throw std::logic_error("Error: " + read.name + " is not cached");
    }
    const ReadType read_type = it->second;
    return read_type;
}

void ReadCache::eraseRead(const Read& read)
{
    const ReadType read_type = typeOfRead(read);

    const auto it = readTypes_.find(read.name);
    readTypes_.erase(it);

    if (read_type == ReadType::kIrrRead || read_type == ReadType::kAnchorRead)
    {
        irrAndAnchorLocations_.erase(read.name);
    }
    if (read_type == ReadType::kIrrRead)
    {
        irrUnits_.erase(read.name);
    }
}

RegionWithCount ReadCache::extractRegionOfIrrOrAnchor(const Read& read) { return irrAndAnchorLocations_.at(read.name); }
string ReadCache::extractUnitOfIrr(const Read& read) { return irrUnits_.at(read.name); }

void ReadCache::cacheAnchorRead(const Read& read)
{
    readTypes_[read.name] = ReadType::kAnchorRead;

    RegionWithCount read_region = createCountableRegion(read.contigId, read.pos, read.pos + 1);
    irrAndAnchorLocations_.emplace(std::make_pair(read.name, read_region));
}

void ReadCache::cacheInrepeatRead(const Read& read, const string& unit)
{
    readTypes_[read.name] = ReadType::kIrrRead;

    RegionWithCount read_region = createCountableRegion(read.contigId, read.pos, read.pos + 1);
    irrAndAnchorLocations_.emplace(std::make_pair(read.name, read_region));

    assert(!unit.empty());
    irrUnits_[read.name] = unit;
}

void ReadCache::cacheOtherRead(const Read& read) { readTypes_[read.name] = ReadType::kOtherRead; }

string ReadCache::printStats()
{
    const string stats = "Cache stats: # reads = " + to_string(readTypes_.size()) + "; # irr and anchor regions = "
        + to_string(irrAndAnchorLocations_.size()) + "; # repeat units = " + to_string(irrUnits_.size());
    return stats;
}

void PairCollector::addAnchor(const Read& read)
{
    if (unparedCache_.isReadCached(read))
    {
        const ReadType mate_type = unparedCache_.typeOfRead(read);
        if (mate_type == ReadType::kIrrRead)
        {
            RegionWithCount irr_region = unparedCache_.extractRegionOfIrrOrAnchor(read);
            const string irr_unit = unparedCache_.extractUnitOfIrr(read);

            RegionWithCount anchor_region = createCountableRegion(read.contigId, read.pos, read.pos + 1);
            anchorRegions_[irr_unit].push_back(anchor_region);
            irrRegions_[irr_unit].push_back(irr_region);
        }
        unparedCache_.eraseRead(read);
    }
    else
    {
        unparedCache_.cacheAnchorRead(read);
    }
}

void PairCollector::addIrr(const Read& read, const std::string& unit)
{
    if (unparedCache_.isReadCached(read))
    {
        const ReadType mate_type = unparedCache_.typeOfRead(read);
        if (mate_type == ReadType::kIrrRead)
        {
            RegionWithCount mate_region = unparedCache_.extractRegionOfIrrOrAnchor(read);
            const string mate_unit = unparedCache_.extractUnitOfIrr(read);

            if (unit == mate_unit)
            {
                RegionWithCount read_region = createCountableRegion(read.contigId, read.pos, read.pos + 1);
                irrRegions_[unit].push_back(read_region);
                irrRegions_[unit].push_back(mate_region);
            }
        }
        else if (mate_type == ReadType::kAnchorRead)
        {
            RegionWithCount irr_region = createCountableRegion(read.contigId, read.pos, read.pos + 1);
            irrRegions_[unit].push_back(irr_region);

            RegionWithCount mate_region = unparedCache_.extractRegionOfIrrOrAnchor(read);
            anchorRegions_[unit].push_back(mate_region);
        }
        unparedCache_.eraseRead(read);
    }
    else
    {
        unparedCache_.cacheInrepeatRead(read, unit);
    }
}

void PairCollector::addOtherRead(const Read& read)
{
    if (unparedCache_.isReadCached(read))
    {
        unparedCache_.eraseRead(read);
    }
    else
    {
        unparedCache_.cacheOtherRead(read);
    }
}

string PairCollector::PrintStats()
{
    string stats = "Collector stats: # anchor regions = " + to_string(anchorRegions_.size()) + "; # irr regions "
        + to_string(irrRegions_.size());

    stats += " " + unparedCache_.printStats();
    return stats;
}
