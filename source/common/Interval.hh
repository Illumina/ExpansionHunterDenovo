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

#include <stdexcept>
#include <string>

class Interval
{
public:
    Interval(int start, int end)
        : start_(start)
        , end_(end)
    {
        if (start_ > end_)
        {
            const auto interval = "(" + std::to_string(start_) + ", " + std::to_string(end_) + ")";
            throw std::runtime_error("Invalid interval endpoints " + interval);
        }
    }

    int start() const { return start_; }
    int end() const { return end_; }
    bool contains(int value) const { return start_ <= value && value <= end_; }

private:
    int start_;
    int end_;
};
