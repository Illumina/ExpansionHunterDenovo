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
