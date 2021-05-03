//
// ExpansionHunter Denovo
// Copyright 2016-2021 Illumina, Inc.
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

#import "common/Manifest.hh"

#include <fstream>
#include <sstream>

using std::string;
using std::vector;

SampleStatus decodeSampleStatus(const string& encoding)
{
    if (encoding == "case")
    {
        return SampleStatus::kCase;
    }
    else if (encoding == "control")
    {
        return SampleStatus::kControl;
    }
    else
    {
        throw std::runtime_error(encoding + " is not a valid sample status");
    }
}

ManifestEntry::ManifestEntry(SampleStatus status, string path)
    : status(status)
    , path(std::move(path))
{
}

ManifestEntry::ManifestEntry(const string& statusEncoding, string path)
    : status(decodeSampleStatus(statusEncoding))
    , path(std::move(path))
{
}

Manifest loadManifest(const string& path, vector<string>& orderedSamples)
{
    Manifest manifest;
    orderedSamples.clear();

    std::ifstream manifestFile(path);
    if (!manifestFile)
    {
        throw std::runtime_error("Unable to load manifest from " + path);
    }

    string line;
    while (std::getline(manifestFile, line))
    {
        std::istringstream decoder(line);
        string sample;
        string statusEncoding;
        string path;
        if (!(decoder >> sample >> statusEncoding >> path))
        {
            throw std::runtime_error("Unable to decode manifest line " + line);
        }
        orderedSamples.push_back(sample);
        manifest.emplace(sample, ManifestEntry(statusEncoding, path));
    }

    return manifest;
}