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

#include "common/Interval.hh"
#include "profile/PairCollector.hh"

ReadType
classifyRead(Interval motifSizeRange, int max_irr_mapq, int min_anchor_mapq, const Read& read, std::string& unit);
PairType
classifyPair(ReadType read_type, const std::string& read_unit, ReadType mate_type, const std::string& mate_unit);
