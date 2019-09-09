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

#include <string>
#include <vector>

class MergeWorkflowParameters
{
public:
    MergeWorkflowParameters(
        const std::string& pathToReference, const std::string& outputPrefix, std::string pathToManifest,
        int shortestUnitToConsider, int longestUnitToConsider);

    const std::string& pathToReference() const { return pathToReference_; }
    const std::string& pathToMultisampleProfile() const { return pathToMultisampleProfile_; }
    const std::string& pathToManifest() const { return pathToManifest_; }
    int shortestUnitToConsider() const { return shortestUnitToConsider_; }
    int longestUnitToConsider() const { return longestUnitToConsider_; }

private:
    std::string pathToReference_;
    std::string pathToMultisampleProfile_;
    std::string pathToManifest_;
    int shortestUnitToConsider_;
    int longestUnitToConsider_;
};

void assertValidity(const MergeWorkflowParameters& parameters);
