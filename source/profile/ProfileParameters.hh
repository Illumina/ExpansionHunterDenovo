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

#include <string>

#include <boost/optional.hpp>

#include "common/Interval.hh"

class ProfileWorkflowParameters
{
public:
    ProfileWorkflowParameters(
        const std::string& outputPrefix, bool logReads, std::string pathToReads, std::string pathToReference,
        Interval motifSizeRange, int minMapqOfAnchorRead, int maxMapqOfInrepeatRead);

    const std::string& profilePath() const { return profilePath_; }
    const std::string& pathToLocusTable() const { return pathToLocusTable_; }
    const std::string& pathToMotifTable() const { return pathToMotifTable_; }
    const std::string& pathToReads() const { return pathToReads_; }
    const std::string& pathToReference() const { return pathToReference_; }
    const boost::optional<std::string>& pathToReadLog() const { return pathToReadLog_; }
    const Interval& motifSizeRange() const { return motifSizeRange_; }
    int minMapqOfAnchorRead() const { return minMapqOfAnchorRead_; }
    int maxMapqOfInrepeatRead() const { return maxMapqOfInrepeatRead_; }

private:
    std::string profilePath_;
    std::string pathToLocusTable_;
    std::string pathToMotifTable_;
    std::string pathToReads_;
    std::string pathToReference_;
    boost::optional<std::string> pathToReadLog_;
    Interval motifSizeRange_;
    int minMapqOfAnchorRead_;
    int maxMapqOfInrepeatRead_;
};

void assertValidity(const ProfileWorkflowParameters& parameters);
