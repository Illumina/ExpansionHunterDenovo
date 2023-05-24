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

#include "ReadClassification.hh"
#include "reads/IrrFinder.hh"

using std::string;

ReadType classifyRead(Interval motifSizeRange, int max_irr_mapq, int min_anchor_mapq, const Read& read, string& unit)
{
    const bool is_unmapped = read.flag & 0x4;
    const bool is_low_mapq = read.mapq <= max_irr_mapq;

    const bool is_irr = (is_unmapped || is_low_mapq) && IsInrepeatRead(read.bases, read.quals, unit, motifSizeRange);

    if (is_irr)
    {
        return ReadType::kIrrRead;
    }

    if (read.mapq >= min_anchor_mapq)
    {
        return ReadType::kAnchorRead;
    }

    return ReadType::kOtherRead;
}

PairType classifyPair(ReadType read_type, const string& read_unit, ReadType mate_type, const string& mate_unit)
{
    if ((read_type == ReadType::kAnchorRead && mate_type == ReadType::kIrrRead)
        || (read_type == ReadType::kIrrRead && mate_type == ReadType::kAnchorRead))
    {
        return PairType::kIrrAnchorPair;
    }

    if (read_type == ReadType::kIrrRead && mate_type == ReadType::kIrrRead && read_unit == mate_unit)
    {
        return PairType::kIrrIrrPair;
    }

    return PairType::kOtherPair;
}
