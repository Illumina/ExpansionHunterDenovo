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

// Code to detect sequences (in-repeat reads) that consist of repetitions of some short string (repeat unit). For
// example, CGGCGGCGGCGG is an in-repeat read with repeat unit CGG.

#pragma once

#include <string>

#include "common/Interval.hh"

int MaxMatchesAtOffset(int offset, const std::string& bases);
double MatchFrequencyAtOffset(int offset, const std::string& bases);
int SmallestFrequentPeriod(
    double minFrequency, const std::string& bases, const Interval& periodSizeRange = Interval(1, 20));
char ExtractConsensusBase(int32_t offset, int32_t period, const std::string& bases);
std::string ExtractConsensusRepeatUnit(double period, const std::string& bases);
std::string MinimialUnitUnderShift(const std::string& unit);
std::string ComputeCanonicalRepeatUnit(const std::string& unit);
std::string ComputeCanonicalRepeatUnit(
    double minFrequency, const std::string& bases, const Interval& motifSizeRange = Interval(1, 20));
bool IsInrepeatRead(
    const std::string& bases, const std::string& quals, std::string& unit,
    const Interval& motifSizeRange = Interval(1, 20));