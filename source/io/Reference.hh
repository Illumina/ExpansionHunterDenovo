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

#include "htslib/faidx.h"

#include "region/GenomicRegion.hh"
#include "region/ReferenceContigInfo.hh"

class Reference
{
public:
    explicit Reference(const std::string& referencePath);
    ~Reference();

    std::string getSequence(const std::string& contigIndex, int64_t start, int64_t end) const;
    std::string getSequence(const GenomicRegion& region) const;

    const ReferenceContigInfo& contigInfo() const { return contigInfo_; }

private:
    std::string referencePath_;
    faidx_t* htsFastaIndexPtr_;
    ReferenceContigInfo contigInfo_;
};
