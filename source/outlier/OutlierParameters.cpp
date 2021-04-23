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

#include "outlier/OutlierParameters.hh"

using boost::optional;
using std::string;

OutlierWorkflowParameters::OutlierWorkflowParameters(
    const string& outputPrefix, string pathToManifest, string pathToMultisampleProfile,
    optional<string> pathToTargetRegions)
    : pathToManifest_(std::move(pathToManifest))
    , pathToMultisampleProfile_(std::move(pathToMultisampleProfile))
    , pathToTargetRegions_(std::move(pathToTargetRegions))
    , pathToLocusAnalysis_(outputPrefix + ".outlier_locus.tsv")
    , pathToMotifAnalysis_(outputPrefix + ".outlier_motif.tsv")
{
}
