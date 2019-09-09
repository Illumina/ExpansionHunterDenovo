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

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "region/ReferenceContigInfo.hh"

class SampleRunStats
{
public:
    SampleRunStats(int meanReadLength, double depth)
        : meanReadLength_(meanReadLength)
        , depth_(depth)
    {
    }

    int meanReadLength() const { return meanReadLength_; }
    double depth() const { return depth_; }

    bool operator==(const SampleRunStats& other) const;

private:
    int meanReadLength_;
    double depth_;
};

std::ostream& operator<<(std::ostream& out, const SampleRunStats& stats);

// Computes read and coverage statistics for each locus from reads aligning to the flanks
class SampleRunStatsCalculator
{
public:
    explicit SampleRunStatsCalculator(ReferenceContigInfo contigInfo);

    void inspect(int contigId, int readLength);

    boost::optional<SampleRunStats> estimate() const;

private:
    ReferenceContigInfo contigInfo_;

    std::unordered_map<int, int64_t> contigIdToReadCount;
    int64_t totalReadCount;
    int64_t sumOfReadLengths;
};
