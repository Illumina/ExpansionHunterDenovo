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
#include <vector>

std::vector<std::vector<std::string>> ShiftUnits(const std::vector<std::string>& units);

double MatchRepeatRc(
    const std::vector<std::vector<std::string>>& units_shifts, const std::string& bases, const std::string& quals,
    size_t min_baseq = 20);

double MatchRepeat(
    const std::vector<std::vector<std::string>>& units_shifts, const std::string& bases, const std::string& quals,
    size_t& match_offset, size_t min_baseq = 20);

double MatchRepeat(
    const std::vector<std::string>& units, const std::string& bases, const std::string& quals, size_t min_baseq = 20);

double MatchUnits(
    const std::vector<std::string>& units, std::string::const_iterator bases_start,
    std::string::const_iterator bases_end, std::string::const_iterator quals_start,
    std::string::const_iterator quals_end, size_t min_baseq = 20);

std::vector<std::vector<std::string>> ShiftUnits(const std::vector<std::string>& units);
