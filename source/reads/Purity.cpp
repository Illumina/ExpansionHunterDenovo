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
#include "reads/Purity.hh"

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "common/SequenceUtils.hh"

using std::string;
using std::vector;

double
MatchRepeatRc(const vector<vector<string>>& units_shifts, const string& bases, const string& quals, size_t min_baseq)
{
    size_t best_offset_forward = 0;
    double forward_score = MatchRepeat(units_shifts, bases, quals, best_offset_forward, min_baseq);
    const string bases_rc = reverseComplement(bases);
    string quals_rc = quals;
    std::reverse(quals_rc.begin(), quals_rc.end());
    size_t best_offset_reverse = 0;
    double reverse_score = MatchRepeat(units_shifts, bases_rc, quals_rc, best_offset_reverse, min_baseq);

    return std::max(forward_score, reverse_score);
}

double MatchRepeat(
    const vector<vector<string>>& units_shifts, const string& bases, const string& quals, size_t& match_offset,
    size_t min_baseq)
{
    double max_score = std::numeric_limits<double>::lowest();
    size_t offset = 0;
    for (const vector<string>& units_shift : units_shifts)
    {
        const double score = MatchRepeat(units_shift, bases, quals, min_baseq);
        if (score > max_score)
        {
            max_score = score;
            match_offset = offset;
        }
        ++offset;
    }
    return max_score;
}

vector<vector<string>> ShiftUnits(const vector<string>& units)
{
    size_t unit_length = units[0].length();
    vector<string> extended_units;
    for (const string& unit : units)
    {
        extended_units.push_back(unit + unit);
    }

    vector<vector<string>> shifted_units;
    for (size_t offset = 0; offset != unit_length; ++offset)
    {
        vector<string> current_shift_units;
        for (const string& extended_unit : extended_units)
        {
            current_shift_units.push_back(extended_unit.substr(offset, unit_length));
        }
        shifted_units.push_back(current_shift_units);
    }

    return shifted_units;
}

double MatchRepeat(const vector<string>& units, const string& bases, const string& quals, size_t min_baseq)
{
    const size_t unit_len = units[0].length();
    double score = 0;
    size_t pos = 0;
    while (pos + unit_len <= bases.length())
    {
        score += MatchUnits(
            units, bases.begin() + pos, bases.begin() + pos + unit_len, quals.begin() + pos,
            quals.begin() + pos + unit_len, min_baseq);
        pos += unit_len;
    }

    if (pos != bases.length())
    {
        score += MatchUnits(units, bases.begin() + pos, bases.end(), quals.begin() + pos, quals.end(), min_baseq);
    }

    return score;
}

double MatchUnits(
    const vector<string>& units, string::const_iterator bases_first, string::const_iterator bases_last,
    string::const_iterator quals_first, string::const_iterator quals_last, size_t min_baseq)
{
    double max_match_count = std::numeric_limits<double>::lowest();
    const size_t kBaseQualOffset = 33;
    const double kMatchScore = 1.0;
    const double kLowqualMismatchScore = 0.5;
    const double kMismatchPenalty = -1.0;

    for (const string& unit : units)
    {
        string::const_iterator unit_first = unit.begin();
        string::const_iterator bases_first_copy = bases_first;
        string::const_iterator quals_first_copy = quals_first;

        double match_count = 0;
        while (bases_first_copy != bases_last)
        {
            if (*bases_first_copy == *unit_first)
            {
                match_count += kMatchScore;
            }
            else if (*quals_first_copy - kBaseQualOffset < min_baseq)
            {
                match_count += kLowqualMismatchScore;
            }
            else
            {
                match_count += kMismatchPenalty;
            }

            ++bases_first_copy;
            ++quals_first_copy;
            ++unit_first;
        }

        if (match_count > max_match_count)
        {
            max_match_count = match_count;
        }
    }

    return max_match_count;
}
