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

#include "classification/irr_finder.h"

#include <string>
#include <vector>

#include "gmock/gmock.h"

using std::string;
using std::vector;
using testing::DoubleNear;

TEST(MaxMatches, CalculatedForVariousOffsets)
{
    const string bases = "ATCGATCG";
    const size_t num_bases = bases.length();
    EXPECT_EQ(num_bases, MaxMatchesAtOffset(0, bases));
    EXPECT_EQ(num_bases - 1, MaxMatchesAtOffset(1, bases));
    EXPECT_EQ(num_bases - 2, MaxMatchesAtOffset(2, bases));
    EXPECT_EQ(0, MaxMatchesAtOffset(num_bases, bases));
    EXPECT_EQ(0, MaxMatchesAtOffset(num_bases + 1, bases));
}

TEST(MatchFrequency, CalculatedForValidOffsets)
{
    const string bases = "GGCCCCGGCCCC";
    vector<double> expected_frequencies = { 0.73, 0.40, 0.33, 0.25, 0.57, 1.00 };

    for (int32_t offset = 1; offset != 7; ++offset)
    {
        const double expected_frequency = expected_frequencies[offset - 1];
        EXPECT_THAT(MatchFrequencyAtOffset(offset, bases), DoubleNear(expected_frequency, 0.005));
    }
}

TEST(MatchFrequency, CalculatedForImperfectRepeat)
{
    const string bases = "ATGATCATGTTGATG";
    const double expected_frequency = 8 / 12.0;
    EXPECT_THAT(MatchFrequencyAtOffset(3, bases), DoubleNear(expected_frequency, 0.005));
}

TEST(MatchFrequency, ExceptionThrownForInvalidOffsets)
{
    const string bases = "GGCCCCGGCCCC";
    const int32_t num_bases = bases.length();
    EXPECT_ANY_THROW(MatchFrequencyAtOffset(-1, bases));
    EXPECT_ANY_THROW(MatchFrequencyAtOffset(0, bases));
    EXPECT_ANY_THROW(MatchFrequencyAtOffset(num_bases / 2 + 1, bases));
}

TEST(SmallestFrequentPeriod, CalculatedForTypicalSequences)
{
    {
        const double min_frequency = 0.85;
        const string bases = "GGCCCCGGCCCC";
        const int32_t period = SmallestFrequentPeriod(min_frequency, bases);
        EXPECT_EQ(6, period);
    }

    {
        const double min_frequency = 0.85;
        const string bases = "ATGATCATGATGATGATGATG";
        const int32_t period = SmallestFrequentPeriod(min_frequency, bases);
        EXPECT_EQ(3, period);
    }
}

TEST(SmallestFrequentPeriod, ReturnsNegativeOneIfNoFrequentPeriodExists)
{
    const double min_frequency = 0.85;
    const string bases = "ATCGGCTA";
    const int32_t period = SmallestFrequentPeriod(min_frequency, bases);
    EXPECT_EQ(-1, period);
}

TEST(ConsensusBase, ExtractedForTypicalSequence)
{
    const string bases = "CGATGACTG";
    const int32_t period = 3;
    EXPECT_EQ('C', ExtractConsensusBase(0, period, bases));
    EXPECT_EQ('G', ExtractConsensusBase(1, period, bases));
    EXPECT_EQ('A', ExtractConsensusBase(2, period, bases));
}

TEST(ConsensusRepeatUnit, ExtractedForTypicalSequences)
{
    {
        const string bases = "CGGCGGCGG";
        const int32_t period = 3;
        EXPECT_EQ("CGG", ExtractConsensusRepeatUnit(period, bases));
    }
    {
        const string bases = "CGGATTATTATTCGG";
        const int32_t period = 3;
        EXPECT_EQ("ATT", ExtractConsensusRepeatUnit(period, bases));
    }
}

TEST(MinimialUnitUnderShift, ComputedForTypicalUnits)
{
    const string repeat_unit = "GGC";
    EXPECT_EQ("CGG", MinimialUnitUnderShift(repeat_unit));
}

TEST(CanonicalUnit, ComputedForTypicalUnits)
{
    {
        const string repeat_unit = "CGG";
        EXPECT_EQ("CCG", ComputeCanonicalRepeatUnit(repeat_unit));
    }
    {
        const string repeat_unit = "GCC";
        EXPECT_EQ("CCG", ComputeCanonicalRepeatUnit(repeat_unit));
    }
}

TEST(ComputeCanonicalRepeatUnitForRead, TypicalInrepeatReadsDetected)
{
    {
        const string bases = "CGGCGCCGGCGG";
        const double min_frequency = 0.8;
        EXPECT_EQ("CCG", ComputeCanonicalRepeatUnit(min_frequency, bases));
        EXPECT_EQ("", ComputeCanonicalRepeatUnit(min_frequency + 0.05, bases));
    }

    {
        const string bases = "ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA";
        const double min_frequency = 0.8;
        EXPECT_EQ("AACCCC", ComputeCanonicalRepeatUnit(min_frequency, bases));
    }
}

TEST(ComputeCanonicalRepeatUnitForRead, HomopolymerReadsDetected)
{
    const string bases = "CCCCCCC";
    const double min_frequency = 1.0;
    EXPECT_EQ("C", ComputeCanonicalRepeatUnit(min_frequency, bases));
}

TEST(IsInrepeatRead, DeterminedForTypicalReads)
{
    {
        const string irr = "CCCCC";
        const string quals = "$$$$$";
        string unit;
        EXPECT_TRUE(IsInrepeatRead(irr, quals, unit));
        EXPECT_EQ("C", unit);
    }
    {
        const string not_irr = "AAAAACCCCC";
        const string quals = "$$$$$$$$$$";
        string unit;
        EXPECT_FALSE(IsInrepeatRead(not_irr, quals, unit));
    }
    {
        const string irr = "TCCACCCACCTCACCCCCCCCCCCCCCCGCCCCCCCCCCACCCCCCCCGCCCCCCCCCCCGGCCCCCCACTCCCCCCCCCCGGTCCTCCCC"
                           "CCCCCCCACCCTCCCCCCC"
                           "CGCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCC";
        const string quals = "------7----7-----7-777-7-F<--777F777F<J-7--7-7-A7-AFJA<<A-<<-7--7A77---7A-77A77A7---7-7-"
                             "7--77-7-77-777---7<7A<"
                             "A-7A)-)-<)7))77A<JJF))--A<F-)-<-)<---7<J";
        string unit;
        EXPECT_FALSE(IsInrepeatRead(irr, quals, unit));
        EXPECT_EQ("ACCCTCCCCCCCCGCCCCCCCCCCCCCCCCCCTCCCCCCCCCTCCCCCCCCCCGGTCCTCCCCCCCCCCC", unit);
    }
    {
        const string irr = "TCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT"
                           "CATTTCATTTCATTTCATT"
                           "TCATTTCTTTTTTTTTATTTTTTTTTATTTTATATCGGAT";
        const string quals = "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                             "((((((((((((((((((((("
                             "((((((((((((((((((((((((((((((((((((((((";

        string unit;
        EXPECT_TRUE(IsInrepeatRead(irr, quals, unit));
        EXPECT_EQ("AAATG", unit);
    }
}

TEST(IsInreapeatRead, DoesNotConsiderReadsMadeOfNsAsInrepeat)
{
    const string n_bases = "NNNNN";
    const string quals = "$$$$$";
    string unit;
    EXPECT_FALSE(IsInrepeatRead(n_bases, quals, unit));
}
