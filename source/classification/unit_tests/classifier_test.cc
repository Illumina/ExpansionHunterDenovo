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

#include "classification/classifiers.cc"

#include "gtest/gtest.h"

using std::string;

TEST(ReadClassification, IrrThatIsNotAnchor_ClassifiedIrr)
{
    const int max_irr_mapq = 10;
    const int min_anchor_mapq = 50;
    Read read;
    read.bases = "CCCCC";
    read.quals = "$$$$$";
    read.mapq = 0;
    string unit;
    ASSERT_EQ(ReadType::kIrrRead, ClassifyRead(max_irr_mapq, min_anchor_mapq, read, unit));
}

TEST(ReadClassification, IrrThatIsAnchor_ClassifiedIrr)
{
    const int max_irr_mapq = 60;
    const int min_anchor_mapq = 50;
    Read read;
    read.bases = "CCCCC";
    read.quals = "$$$$$";
    read.mapq = 60;
    string unit;
    ASSERT_EQ(ReadType::kIrrRead, ClassifyRead(max_irr_mapq, min_anchor_mapq, read, unit));
}

TEST(ReadClassification, AnchorThatIsNotIrr_ClassifiedAnchor)
{
    const int max_irr_mapq = 10;
    const int min_anchor_mapq = 50;
    Read read;
    read.bases = "ATCGC";
    read.quals = "$$$$$";
    read.mapq = 60;
    string unit;
    ASSERT_EQ(ReadType::kAnchorRead, ClassifyRead(max_irr_mapq, min_anchor_mapq, read, unit));
}

TEST(ReadClassification, NotAnchorAndNotIrr_ClassifiedOther)
{
    const int max_irr_mapq = 10;
    const int min_anchor_mapq = 50;
    Read read;
    read.bases = "ATCGC";
    read.quals = "$$$$$";
    read.mapq = 0;
    string unit;
    ASSERT_EQ(ReadType::kOtherRead, ClassifyRead(max_irr_mapq, min_anchor_mapq, read, unit));
}
