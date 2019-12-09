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

#include "common/Interval.hh"
#include "profile/PairCollector.hh"

ReadType
classifyRead(Interval motifSizeRange, int max_irr_mapq, int min_anchor_mapq, const Read& read, std::string& unit);
PairType
classifyPair(ReadType read_type, const std::string& read_unit, ReadType mate_type, const std::string& mate_unit);