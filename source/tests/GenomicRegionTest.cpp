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

#include "region/GenomicRegion.hh"

#include "thirdparty/catch2/catch.hpp"

#include <limits>

using Catch::Contains;
using std::vector;

TEST_CASE("Region must always be in a valid state", "[region initialization]")
{
    REQUIRE_THROWS_WITH(GenomicRegion(-2, 1, 10), Contains("is invalid"));
}

TEST_CASE("Overlapping regions have zero distance", "[manipulating regions]")
{
    GenomicRegion regionA(1, 1, 10);
    GenomicRegion regionB(1, 5, 15);

    REQUIRE(regionA.distance(regionB) == 0);
}

TEST_CASE("Disjoint regions have non-zero distance", "[manipulating regions]")
{
    GenomicRegion regionA(1, 50, 70);
    GenomicRegion regionB(1, 0, 20);

    REQUIRE(regionA.distance(regionB) == 30);
    REQUIRE(regionB.distance(regionA) == 30);
}

TEST_CASE("Regions on different contigs have maximum possible distance", "[manipulating regions]")
{
    GenomicRegion regionA(1, 50, 70);
    GenomicRegion regionB(2, 0, 20);

    REQUIRE(regionA.distance(regionB) == std::numeric_limits<int64_t>::max());
}

TEST_CASE("Overlapping sorted regions are correctly merged", "[manipulating regions]")
{
    vector<RegionWithCount> regions
        = { createCountableRegion(1, 10, 20), createCountableRegion(1, 15, 25), createCountableRegion(1, 20, 35) };

    sortAndMerge(regions);

    vector<RegionWithCount> expectedRegions = { { 1, 10, 35, CountFeature(3) } };
    REQUIRE(regions == expectedRegions);
}

TEST_CASE("Overlapping unsorted regions are correctly merged", "[manipulating regions]")
{
    vector<RegionWithSampleCount> regions = { { 1, 10, 20, SampleCountFeature({ { "A", 5 }, { "B", 4 } }) },
                                              { 1, 15, 25, SampleCountFeature({ { "A", 3 }, { "B", 4 } }) },
                                              { 1, 20, 35, SampleCountFeature({ { "B", 2 }, { "C", 8 } }) } };

    sortAndMerge(regions);

    vector<RegionWithSampleCount> expectedRegions
        = { RegionWithSampleCount(1, 10, 35, SampleCountFeature({ { "A", 8 }, { "B", 10 }, { "C", 8 } })) };
    REQUIRE(regions == expectedRegions);
}

TEST_CASE("Disjoint regions are correctly merged", "[manipulating regions]")
{
    vector<RegionWithCount> regions
        = { createCountableRegion(1, 15, 25), createCountableRegion(2, 10, 20), createCountableRegion(1, 20, 35) };

    sortAndMerge(regions);

    vector<RegionWithCount> expectedRegions = { { 1, 15, 35, CountFeature(2) }, { 2, 10, 20, CountFeature(1) } };
    REQUIRE(regions == expectedRegions);
}

TEST_CASE("Proximal regions are correctly merged", "[manipulating regions]")
{
    vector<RegionWithCount> regions
        = { createCountableRegion(1, 200, 250), createCountableRegion(1, 500, 550), createCountableRegion(1, 0, 10),
            createCountableRegion(1, 1100, 1200), createCountableRegion(2, 1100, 1200) };

    sortAndMerge(regions);

    vector<RegionWithCount> expectedRegions
        = { { 1, 0, 550, CountFeature(3) }, { 1, 1100, 1200, CountFeature(1) }, { 2, 1100, 1200, CountFeature(1) } };
    REQUIRE(regions == expectedRegions);
}

TEST_CASE("Nested regions are correctly merged", "[manipulating regions]")
{
    vector<RegionWithCount> regions = { createCountableRegion(1, 100, 200), createCountableRegion(1, 90, 300) };

    sortAndMerge(regions);

    vector<RegionWithCount> expectedRegions = { { 1, 90, 300, CountFeature(2) } };
    REQUIRE(regions == expectedRegions);
}

TEST_CASE("Genomic region can be decoded from string", "[manipulating regions]")
{
    ReferenceContigInfo contigInfo({ { "chr1", 0 } });
    GenomicRegion region = decode(contigInfo, "chr1:1-100");

    GenomicRegion expectedRegion(0, 1, 100);
    REQUIRE(region == expectedRegion);
}

TEST_CASE("Genomic regions with ambiguous contigs can be decoded", "[manipulating regions]")
{
    ReferenceContigInfo contigInfo({ { "HLA-DQA1*05:11", 1000 } });
    GenomicRegion region = decode(contigInfo, "HLA-DQA1*05:11:6177-6177");

    GenomicRegion expectedRegion(0, 6177, 6177);
    REQUIRE(region == expectedRegion);
}

TEST_CASE("Unaligned regions can be decoded", "[manipulating regions]")
{
    ReferenceContigInfo contigInfo({ { "chr1", 1000 } });
    GenomicRegion region = decode(contigInfo, "unaligned");

    GenomicRegion expectedRegion(-1, 0, 0);
    REQUIRE(region == expectedRegion);
}