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
