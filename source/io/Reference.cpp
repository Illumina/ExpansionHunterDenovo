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

#include "io/Reference.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

using std::string;
using std::to_string;
using std::vector;

Reference::Reference(const string& referencePath)
    : referencePath_(referencePath)
    , contigInfo_({})
{
    htsFastaIndexPtr_ = fai_load(referencePath_.c_str());

    std::vector<std::pair<std::string, int64_t>> internalNamesAndSizes;

    for (int contigIndex = 0; contigIndex != faidx_nseq(htsFastaIndexPtr_); ++contigIndex)
    {
        const char* sequenceName = faidx_iseq(htsFastaIndexPtr_, contigIndex);
        int64_t sequenceLength = faidx_seq_len(htsFastaIndexPtr_, sequenceName);
        internalNamesAndSizes.emplace_back(sequenceName, sequenceLength);
    }

    contigInfo_ = ReferenceContigInfo(internalNamesAndSizes);
}

Reference::~Reference() { fai_destroy(htsFastaIndexPtr_); }

string Reference::getSequence(const string& contigName, int64_t start, int64_t end) const
{
    const int contigIndex = contigInfo_.getContigId(contigName);
    const char* contigNamePtr = faidx_iseq(htsFastaIndexPtr_, contigIndex);

    int extractedLength;
    // This htslib function is 0-based closed but our coordinates are half open
    char* sequencePtr = faidx_fetch_seq(htsFastaIndexPtr_, contigNamePtr, start, end - 1, &extractedLength);

    if (!sequencePtr || extractedLength < 0 || extractedLength < end - start)
    {
        const string encoding(contigName + ":" + to_string(start) + "-" + to_string(end));
        const string message = "Unable to extract " + encoding + " from " + referencePath_;
        throw std::runtime_error(message);
    }

    string sequence("N", extractedLength);
    std::transform(sequencePtr, sequencePtr + extractedLength, sequence.begin(), ::toupper);
    free(sequencePtr);

    return sequence;
}

string Reference::getSequence(const GenomicRegion& region) const
{
    return getSequence(contigInfo_.getContigName(region.contigId()), region.start(), region.end());
}
