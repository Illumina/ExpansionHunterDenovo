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

#include "reads/Read.hh"

std::ostream& operator<<(std::ostream& out, const Read& read)
{
    out << read.name << "\t" << read.bases << "\t" << read.quals << "\t" << read.contigId << "\t" << read.pos << "\t"
        << read.mateContigId << "\t" << read.matePos << "\t" << read.mapq << "\t" << read.flag;
    return out;
}