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

#include <memory>
#include <string>

class PathParameters
{
public:
    PathParameters(std::string reads, std::string reference, std::string outputPrefix)
        : reads_(std::move(reads))
        , reference_(std::move(reference))
        , outputPrefix_(std::move(outputPrefix))
    {
    }

    const std::string& reads() const { return reads_; };
    const std::string& reference() const { return reference_; };
    const std::string& outputPrefix() const { return outputPrefix_; }

private:
    std::string reads_;
    std::string reference_;
    std::string outputPrefix_;
};

class HeuristicParameters
{
public:
    HeuristicParameters(
        int shortestUnitToConsider, int longestUnitToConsider, int minMapqOfAnchorRead, int maxMapqOfInrepeatRead)
        : shortestUnitToConsider_(shortestUnitToConsider)
        , longestUnitToConsider_(longestUnitToConsider)
        , minMapqOfAnchorRead_(minMapqOfAnchorRead)
        , maxMapqOfInrepeatRead_(maxMapqOfInrepeatRead)
    {
    }

    int shortestUnitToConsider() const { return shortestUnitToConsider_; }
    int longestUnitToConsider() const { return longestUnitToConsider_; }
    int minMapqOfAnchorRead() const { return minMapqOfAnchorRead_; }
    int maxMapqOfInrepeatRead() const { return maxMapqOfInrepeatRead_; }

private:
    int shortestUnitToConsider_;
    int longestUnitToConsider_;
    int minMapqOfAnchorRead_;
    int maxMapqOfInrepeatRead_;
};

class ProgramParameters
{
public:
    ProgramParameters(PathParameters paths, int readLength, HeuristicParameters heuristics)
        : paths_(std::move(paths))
        , readLength_(readLength)
        , heuristics_(heuristics)
    {
    }

    const PathParameters& paths() const { return paths_; }
    int readLength() const { return readLength_; }
    const HeuristicParameters& heuristics() const { return heuristics_; }

private:
    PathParameters paths_;
    int readLength_;
    HeuristicParameters heuristics_;
};