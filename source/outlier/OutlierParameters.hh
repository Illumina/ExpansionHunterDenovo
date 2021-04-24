//
// ExpansionHunter Denovo
// Copyright 2016-2020 Illumina, Inc.
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

#include <boost/optional.hpp>

class OutlierWorkflowParameters
{
public:
    OutlierWorkflowParameters(
        const std::string& outputPrefix, std::string pathToReference, std::string pathToManifest,
        std::string pathToMultisampleProfile, boost::optional<std::string> pathToTargetRegions);

    const std::string& pathToReference() const { return pathToReference_; }
    const std::string& pathToMultisampleProfile() const { return pathToMultisampleProfile_; }
    const std::string& pathToManifest() const { return pathToManifest_; }

private:
    std::string pathToReference_;
    std::string pathToMultisampleProfile_;
    std::string pathToManifest_;
    boost::optional<std::string> pathToTargetRegions_;

    std::string pathToLocusAnalysis_;
    std::string pathToMotifAnalysis_;
};