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
