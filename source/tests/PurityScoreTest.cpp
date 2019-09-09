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

#include "thirdparty/catch2/catch.hpp"

using std::string;
using std::vector;

// TEST_CASE("Reverse complement can be computed for any sequence", "[reverse complementing]")
//{
//    const string bases = "ATCGN";
//    REQUIRE("NCGAT" == reverseComplement(bases));
//}

TEST_CASE("Unit matches to itself", "[calculating purity scores]")
{
    char qual_chars[] = { 40, 40, 40, 40, 40, 40 };
    string quals = qual_chars;
    string bases = "GGCCCC";
    vector<string> units = { "GGCCCC" };

    REQUIRE(MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()) == Approx(6.0));
}

/*
TEST(TestUnitMatching, MatchesMultipleUnits) {
string quals = "PPPPPP";
string bases = "AACTCC";
vector<string> units = {"GGCCCC", "AACTCC"};

EXPECT_DOUBLE_EQ(MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()), 6.0);
}

TEST(TestUnitMatching, MatchesShortSequence) {
string quals = "PPP";
string bases = "AAC";
vector<string> units = {"GGCCCC", "AACTCC"};

EXPECT_DOUBLE_EQ(MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()), 3.0);
}

TEST(TestUnitMatching, MatchesLowqualBases) {
string quals = "(PP(((";
string bases = "AACCGG";
vector<string> units = {"GGCCCC", "AACTCC"};

EXPECT_DOUBLE_EQ(MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()), 4.5);
}

TEST(TestUnitMatching, ScoreCanBeNegative) {
string quals = "PPPPPP";
string bases = "AACCGG";
vector<string> units = {"ATTTTT", "AATTTT"};

EXPECT_DOUBLE_EQ(MatchUnits(units, bases.begin(), bases.end(), quals.begin(), quals.end()), -2.0);
}

TEST(TestRepeatMatching, RepeatMatches) {
string quals = "PPPPPPPP";
string bases = "ACGATGAC";
vector<string> units = {"AAG", "ACG"};

EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 6.0);
}

TEST(TestRepeatMatching, MotifShorterByOne) {
string quals = "PPPPPPPP";
string bases = "ACGATGAC";
vector<string> units = {"AAAATTT", "ACGATGA"};

EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 6.0);
}

TEST(TestRepeatMatching, EmptySequenceScoresZero) {
string quals;
string bases;
vector<string> units = {"AAG", "ACG"};

EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 0);
}

TEST(TestRepeatMatching, SingletonScoresOne) {
string quals = "B";
string bases = "G";
vector<string> units = {"G"};

EXPECT_DOUBLE_EQ(MatchRepeat(units, bases, quals), 1.0);
}

TEST(TestRepeatMatching, MakeShiftedUnits) {
vector<string> units = {"AAG", "ACG"};
vector<vector<string>> shifted_units = {{"AAG", "ACG"}, {"AGA", "CGA"}, {"GAA", "GAC"}};
ASSERT_EQ(ShiftUnits(units), shifted_units);
}

TEST(TestRepeatMatching, RepeatMatchesWithShift) {
vector<string> units = {"AAG", "ACG"};
vector<vector<string>> units_shifts = ShiftUnits(units);
string quals = "PPPPPPPP";
string bases = "CGACGACG";

size_t best_offset = 0;
ASSERT_DOUBLE_EQ(MatchRepeat(units_shifts, bases, quals, best_offset), 8.0);
}

TEST(TestRepeatMatching, CalculatesBestMatchOffset) {
vector<string> units = {"AAG", "ACG"};
vector<vector<string>> units_shifts = ShiftUnits(units);
string quals = "PPPPPPPP";
string bases = "CGACGACG";
size_t offset = 0;
MatchRepeat(units_shifts, bases, quals, offset);
ASSERT_DOUBLE_EQ(offset, 1);
}

TEST(TestRepeatMatching, RepeatMatchesReverseCompliment) {
vector<string> units = {"AAG", "ACG"};
vector<vector<string>> units_shifts = ShiftUnits(units);
string quals = "((PPPPPP";
string bases = "AATCGTCG";
ASSERT_DOUBLE_EQ(MatchRepeatRc(units_shifts, bases, quals), 7.0);
} */