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

#include "io/HtsHelpers.hh"

#include <utility>

using std::pair;
using std::string;
using std::vector;

string decodeQuals(bam1_t* htsAlignPtr)
{
    string quals;
    uint8_t* htsQualsPtr = bam_get_qual(htsAlignPtr);
    const int readLength = htsAlignPtr->core.l_qseq;
    quals.resize(readLength);

    for (int index = 0; index < readLength; ++index)
    {
        quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
    }

    return quals;
}

string decodeBases(bam1_t* htsAlignPtr)
{
    string bases;
    uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
    const int32_t readLength = htsAlignPtr->core.l_qseq;
    bases.resize(readLength);

    for (int32_t index = 0; index < readLength; ++index)
    {
        bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
    }

    return bases;
}

Read decodeHtsRead(bam1_t* htsAlignPtr)
{
    Read read;

    read.bases = decodeBases(htsAlignPtr);
    read.flag = htsAlignPtr->core.flag;
    read.mapq = htsAlignPtr->core.qual;
    read.name = bam_get_qname(htsAlignPtr);
    read.pos = htsAlignPtr->core.pos;
    read.quals = decodeQuals(htsAlignPtr);
    read.contigId = htsAlignPtr->core.tid;
    read.mateContigId = htsAlignPtr->core.mtid;
    read.matePos = htsAlignPtr->core.mpos;

    return read;
}

bool isPrimaryAlignment(bam1_t* htsAlignPtr)
{
    return !((htsAlignPtr->core.flag & BAM_FSECONDARY) || (htsAlignPtr->core.flag & BAM_FSUPPLEMENTARY));
}

ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr)
{
    vector<pair<string, int64_t>> contigNamesAndSizes;
    const int32_t numContigs = htsHeaderPtr->n_targets;
    contigNamesAndSizes.reserve(numContigs);

    for (int32_t contigIndex = 0; contigIndex != numContigs; ++contigIndex)
    {
        const string contig = htsHeaderPtr->target_name[contigIndex];
        int64_t size = htsHeaderPtr->target_len[contigIndex];
        contigNamesAndSizes.push_back(std::make_pair(contig, size));
    }

    return ReferenceContigInfo(contigNamesAndSizes);
}
