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

#include <memory>
#include <string>
#include <vector>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "reads/Read.hh"
#include "region/ReferenceContigInfo.hh"

class HtsFileStreamer
{
public:
    HtsFileStreamer(std::string htsFilePath, std::string referencePath)
        : htsFilePath_(std::move(htsFilePath))
        , referencePath_(std::move(referencePath))
        , contigInfo_({})
    {
        openHtsFile();
        loadHeader();
        prepareForStreamingAlignments();
    }
    ~HtsFileStreamer();

    const ReferenceContigInfo& contigInfo() const { return contigInfo_; }

    bool trySeekingToNextPrimaryAlignment();

    int currentReadContigId() const;
    int currentReadPosition() const;
    int currentReadLength() const;
    int currentMateContigId() const;
    int currentMatePosition() const;

    bool isStreamingAlignedReads() const;

    Read decodeRead() const;

private:
    enum class Status
    {
        kStreamingReads,
        kFinishedStreaming
    };

    void openHtsFile();
    void loadHeader();
    void prepareForStreamingAlignments();

    std::string htsFilePath_;
    std::string referencePath_;
    ReferenceContigInfo contigInfo_;
    Status status_ = Status::kStreamingReads;

    htsFile* htsFilePtr_ = nullptr;
    bam1_t* htsAlignmentPtr_ = nullptr;
    bam_hdr_t* htsHeaderPtr_ = nullptr;
};
