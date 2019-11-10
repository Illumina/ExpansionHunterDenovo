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
