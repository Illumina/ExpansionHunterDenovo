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
