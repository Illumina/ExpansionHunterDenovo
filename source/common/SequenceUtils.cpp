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

#include "SequenceUtils.hh"

using std::string;

string reverseComplement(const string& bases)
{
    string bases_rc = bases;
    string::reverse_iterator bases_rc_iter = bases_rc.rbegin();

    char complemented_base = ' ';
    for (char base : bases)
    {
        switch (base)
        {
        case 'A':
            complemented_base = 'T';
            break;
        case 'C':
            complemented_base = 'G';
            break;
        case 'G':
            complemented_base = 'C';
            break;
        case 'T':
            complemented_base = 'A';
            break;
        default:
            complemented_base = 'N';
        }
        *bases_rc_iter++ = complemented_base;
    }

    return bases_rc;
}
