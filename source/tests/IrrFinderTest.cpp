//
// ExpansionHunter Denovo
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "reads/IrrFinder.hh"

#include <vector>

#include "thirdparty/catch2/catch.hpp"

using Catch::Contains;
using std::string;
using std::vector;

TEST_CASE("Match count calculated for various offsets", "[Determining motif]")
{
    const string bases = "ATCGATCG";
    const size_t num_bases = bases.length();
    REQUIRE(MaxMatchesAtOffset(0, bases) == num_bases);
    REQUIRE(MaxMatchesAtOffset(1, bases) == num_bases - 1);
    REQUIRE(MaxMatchesAtOffset(2, bases) == num_bases - 2);
    REQUIRE(MaxMatchesAtOffset(num_bases, bases) == 0);
    REQUIRE(MaxMatchesAtOffset(num_bases + 1, bases) == 0);
}

TEST_CASE("Match frequency calculated for various offsets", "[Determining motif]")
{
    const string bases = "GGCCCCGGCCCC";
    vector<double> expected_frequencies = { 0.73, 0.40, 0.33, 0.25, 0.57, 1.00 };

    for (int offset = 1; offset != 7; ++offset)
    {
        const double expected_frequency = expected_frequencies[offset - 1];
        REQUIRE(MatchFrequencyAtOffset(offset, bases) == Approx(expected_frequency).margin(0.01));
    }
}

TEST_CASE("Match frequency calculated for imperfect repeat", "[Determining motif]")
{
    const string bases = "ATGATCATGTTGATG";
    const double expected_frequency = 8 / 12.0;
    REQUIRE(MatchFrequencyAtOffset(3, bases) == Approx(expected_frequency));
}

TEST_CASE("Match frequency calculation throws for on invalid offsets", "[Determining motif]")
{
    const string bases = "GGCCCCGGCCCC";
    const int num_bases = bases.length();
    REQUIRE_THROWS_WITH(MatchFrequencyAtOffset(-1, bases), Contains("not a valid offset"));
    REQUIRE_THROWS_WITH(MatchFrequencyAtOffset(0, bases), Contains("not a valid offset"));
    // REQUIRE_THROWS_WITH(MatchFrequencyAtOffset(num_bases / 2 + 1, bases), Contains("not a valid offset"));
}

TEST_CASE("Calculating period for typical sequences", "[Determining motif]")
{
    {
        const double min_frequency = 0.85;
        const string bases = "GGCCCCGGCCCC";
        const int period = SmallestFrequentPeriod(min_frequency, bases);
        REQUIRE(period == 6);
    }

    {
        const double min_frequency = 0.85;
        const string bases = "ATGATCATGATGATGATGATG";
        const int period = SmallestFrequentPeriod(min_frequency, bases);
        REQUIRE(period == 6); // 3 might be a better answers
    }
}

TEST_CASE("-1 is used as a special value for no period", "[Determining motif]")
{
    const double min_frequency = 0.85;
    const string bases = "ATCGGCTA";
    const int period = SmallestFrequentPeriod(min_frequency, bases);
    REQUIRE(period == -1);
}

TEST_CASE("Consensus base is extracted a sequence", "[Determining motif]")
{
    const string bases = "CGATGACTG";
    const int period = 3;
    REQUIRE(ExtractConsensusBase(0, period, bases) == 'C');
    REQUIRE(ExtractConsensusBase(1, period, bases) == 'G');
    REQUIRE(ExtractConsensusBase(2, period, bases) == 'A');
}

TEST_CASE("Consensus motif is extracted for a sequence", "[Determining motif]")
{
    {
        const string bases = "CGGCGGCGG";
        const int period = 3;
        REQUIRE(ExtractConsensusRepeatUnit(period, bases) == "CGG");
    }
    {
        const string bases = "CGGATTATTATTCGG";
        const int period = 3;
        REQUIRE(ExtractConsensusRepeatUnit(period, bases) == "ATT");
    }
}

TEST_CASE("Minimal motif under shift is computed a motif", "[Determining motif]")
{
    const string repeat_unit = "GGC";
    REQUIRE(MinimialUnitUnderShift(repeat_unit) == "CGG");
}

TEST_CASE("Canonical motif computed for a motif", "[Determining motif]")
{
    {
        const string repeat_unit = "CGG";
        REQUIRE(ComputeCanonicalRepeatUnit(repeat_unit) == "CCG");
    }
    {
        const string repeat_unit = "GCC";
        REQUIRE(ComputeCanonicalRepeatUnit(repeat_unit) == "CCG");
    }
}

TEST_CASE("ComputeCanonicalRepeatUnitForRead, TypicalInrepeatReadsDetected", "[Determining motif]")
{
    {
        const string bases = "CGGCGCCGGCGG";
        const double min_frequency = 0.8;
        REQUIRE(ComputeCanonicalRepeatUnit(min_frequency, bases) == "CCG");
        REQUIRE(ComputeCanonicalRepeatUnit(min_frequency + 0.05, bases).empty());
    }

    {
        const string bases = "ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA";
        const double min_frequency = 0.8;
        REQUIRE(ComputeCanonicalRepeatUnit(min_frequency, bases) == "AACCCC");
    }
}

TEST_CASE("Compute canonical motif for homopolymer run", "[Determining motif]")
{
    const string bases = "CCCCCCC";
    const double min_frequency = 1.0;
    REQUIRE(ComputeCanonicalRepeatUnit(min_frequency, bases) == "C");
}

TEST_CASE("Irr check is performed on various reads", "[Determining motif]")
{
    {
        const string irr = "CCCCC";
        const string quals = "$$$$$";
        string unit;
        REQUIRE(IsInrepeatRead(irr, quals, unit));
        REQUIRE(unit == "C");
    }
    {
        const string not_irr = "AAAAACCCCC";
        const string quals = "$$$$$$$$$$";
        string unit;
        REQUIRE(!IsInrepeatRead(not_irr, quals, unit));
    }
    {
        const string irr = "TCCACCCACCTCACCCCCCCCCCCCCCCGCCCCCCCCCCACCCCCCCCGCCCCCCCCCCCGGCCCCCCACTCCCCCCCCCCGGTCCTCCCC"
                           "CCCCCCCACCCTCCCCCCCCGCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCC";
        const string quals = "------7----7-----7-777-7-F<--777F777F<J-7--7-7-A7-AFJA<<A-<<-7--7A77---7A-77A77A7---7-7-"
                             "7--77-7-77-777---7<7A<A-7A)-)-<)7))77A<JJF))--A<F-)-<-)<---7<J";
        string unit;
        REQUIRE(!IsInrepeatRead(irr, quals, unit, Interval(1, 70)));
        // REQUIRE(unit == "ACCCTCCCCCCCCGCCCCCCCCCCCCCCCCCCTCCCCCCCCCTCCCCCCCCCCGGTCCTCCCCCCCCCCC");
    }
    {
        const string irr = "TCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT"
                           "CATTTCATTTCATTTCATTTCATTTCTTTTTTTTTATTTTTTTTTATTTTATATCGGAT";
        const string quals = "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                             "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((";

        string unit;
        REQUIRE(IsInrepeatRead(irr, quals, unit));
        REQUIRE(unit == "AAATG");
    }
    {
        string irr = "CCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCCCGCCCCGCC"
                     "CCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCGCCCCGCCCCGCG";
        string quals = "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                       "(((((((((((((((((((((((((((((((((((((((((((((((((((((((";
        string unit;
        Interval motifSizeRange(1, 15);
        REQUIRE(IsInrepeatRead(irr, quals, unit, motifSizeRange));
        REQUIRE(unit == "CCCCGCCCCGCG");
    }
    {
        string irr = "GGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCG"
                     "GGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCG";
        string quals = "((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
                       "(((((((((((((((((((((((((((((((((((((((((((((((((((((((";
        string unit;
        Interval motifSizeRange(1, 20);
        REQUIRE(IsInrepeatRead(irr, quals, unit, motifSizeRange));
        REQUIRE(unit == "CCCCGCCCCGCG");
    }
}

TEST_CASE("Reads made Ns are not inrepeat reads", "[Determining motif]")
{
    const string n_bases = "NNNNN";
    const string quals = "$$$$$";
    string unit;
    REQUIRE(!IsInrepeatRead(n_bases, quals, unit));
}
