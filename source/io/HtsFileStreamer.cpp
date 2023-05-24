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

#include "io/HtsFileStreamer.hh"

#include <stdexcept>

//#include "boost/filesystem.hpp"

#include "io/HtsHelpers.hh"

using std::string;

void HtsFileStreamer::openHtsFile()
{
    htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

    if (!htsFilePtr_)
    {
        throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
    }

    // Set reference index
    const string referenceIndex = referencePath_ + ".fai";
    // if (!boost::filesystem::exists(referenceIndex))
    //{
    //    throw std::runtime_error("Reference index does not exist: " + referenceIndex);
    //}

    if (hts_set_fai_filename(htsFilePtr_, referenceIndex.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set reference index");
    }
    // End setting reference index
}

void HtsFileStreamer::loadHeader()
{
    htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

    if (!htsHeaderPtr_)
    {
        throw std::runtime_error("Failed to read header of " + htsFilePath_);
    }

    contigInfo_ = decodeContigInfo(htsHeaderPtr_);
}

void HtsFileStreamer::prepareForStreamingAlignments() { htsAlignmentPtr_ = bam_init1(); }

bool HtsFileStreamer::trySeekingToNextPrimaryAlignment()
{
    if (status_ != Status::kStreamingReads)
    {
        return false;
    }

    int32_t returnCode = 0;

    while ((returnCode = sam_read1(htsFilePtr_, htsHeaderPtr_, htsAlignmentPtr_)) >= 0)
    {
        if (isPrimaryAlignment(htsAlignmentPtr_))
            return true;
    }

    status_ = Status::kFinishedStreaming;

    if (returnCode < -1)
    {
        throw std::runtime_error("Failed to extract a record from " + htsFilePath_);
    }

    return false;
}

int HtsFileStreamer::currentReadContigId() const { return htsAlignmentPtr_->core.tid; }
int HtsFileStreamer::currentReadPosition() const { return htsAlignmentPtr_->core.pos; }
int HtsFileStreamer::currentMateContigId() const { return htsAlignmentPtr_->core.mtid; }
int HtsFileStreamer::currentMatePosition() const { return htsAlignmentPtr_->core.mpos; }
int HtsFileStreamer::currentReadLength() const { return htsAlignmentPtr_->core.l_qseq; }

bool HtsFileStreamer::isStreamingAlignedReads() const
{
    return status_ != Status::kFinishedStreaming && currentReadContigId() != -1;
}

Read HtsFileStreamer::decodeRead() const { return decodeHtsRead(htsAlignmentPtr_); }

HtsFileStreamer::~HtsFileStreamer()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;
}
