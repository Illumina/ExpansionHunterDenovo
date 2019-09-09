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

#include "merge/MergeParameters.hh"

#include <memory>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

using std::string;

MergeWorkflowParameters::MergeWorkflowParameters(
    const std::string& pathToReference, const string& outputPrefix, string pathToManifest, int shortestUnitToConsider,
    int longestUnitToConsider)
    : pathToReference_(pathToReference)
    , pathToMultisampleProfile_(outputPrefix + ".multisample_profile.json")
    , pathToManifest_(std::move(pathToManifest))
    , shortestUnitToConsider_(shortestUnitToConsider)
    , longestUnitToConsider_(longestUnitToConsider)
{
}

/*

static void assertIndexExists(const string& htsFilePath)
{
    const vector<string> kPossibleIndexExtensions = { ".bai", ".csi", ".crai" };

    for (const string& indexExtension : kPossibleIndexExtensions)
    {
        if (fs::exists(htsFilePath + indexExtension))
        {
            return;
        }
    }

    throw std::invalid_argument("Could not find index of " + htsFilePath);
}

void assertValidity(const UserParameters& userParameters)
{
    // Validate input file paths
    assertPathToExistingFile(userParameters.htsFilePath);
    assertIndexExists(userParameters.htsFilePath);
    assertPathToExistingFile(userParameters.referencePath);
    assertPathToExistingFile(userParameters.catalogPath);

    // Validate output prefix
    assertWritablePath(userParameters.outputPrefix);

    // Validate sample parameters
    if (userParameters.sampleSexEncoding != "female" && userParameters.sampleSexEncoding != "male")
    {
        throw std::invalid_argument(userParameters.sampleSexEncoding + " is not a valid sex encoding");
    }

    // Heuristic parameters
    if (userParameters.alignerType != "dag-aligner" && userParameters.alignerType != "path-aligner")
    {
        throw std::invalid_argument(userParameters.alignerType + " is not a valid aligner type");
    }

    const int kMinExtensionLength = 500;
    const int kMaxExtensionLength = 1500;
    if (userParameters.regionExtensionLength < kMinExtensionLength
        && userParameters.regionExtensionLength > kMaxExtensionLength)
    {
        const string message = "Extension length of size " + to_string(userParameters.regionExtensionLength)
            + " is not supported; the range of allowed extensions is between " + to_string(kMinExtensionLength)
            + " and " + to_string(kMaxExtensionLength);
        throw std::invalid_argument(message);
    }

    const int kMinQualityCutoffForGoodBaseCall = 5;
    const int kMaxQualityCutoffForGoodBaseCall = 40;
    if (userParameters.qualityCutoffForGoodBaseCall < kMinQualityCutoffForGoodBaseCall
        && userParameters.qualityCutoffForGoodBaseCall > kMaxQualityCutoffForGoodBaseCall)
    {
        const string message = "Base call quality cutoff of " + to_string(userParameters.qualityCutoffForGoodBaseCall)
            + " is not supported; the range of allowed cutoffs is between "
            + to_string(kMinQualityCutoffForGoodBaseCall) + " and " + to_string(kMaxQualityCutoffForGoodBaseCall);
        throw std::invalid_argument(message);
    }
}

SampleParameters decodeSampleParameters(const UserParameters& userParams)
{
    fs::path boostHtsFilePath(userParams.htsFilePath);
    auto sampleId = boostHtsFilePath.stem().string();
    Sex sex = decodeSampleSex(userParams.sampleSexEncoding);
    return SampleParameters(sampleId, sex);
}

AnalysisMode decodeAnalysisMode(const string& encoding)
{
    if (encoding == "streaming")
    {
        return AnalysisMode::kStreaming;
    }
    else if (encoding == "seeking")
    {
        return AnalysisMode::kSeeking;
    }
    else
    {
        throw std::logic_error("Invalid encoding of data input mode '" + encoding + "'");
    }
}

LogLevel decodeLogLevel(const string& encoding)
{
    if (encoding == "trace")
    {
        return LogLevel::kTrace;
    }
    if (encoding == "debug")
    {
        return LogLevel::kDebug;
    }
    else if (encoding == "info")
    {
        return LogLevel::kInfo;
    }
    else if (encoding == "warn")
    {
        return LogLevel::kWarn;
    }
    else if (encoding == "error")
    {
        return LogLevel::kError;
    }
    else
    {
        throw std::logic_error("Invalid encoding of logging level " + encoding);
    }
}

boost::optional<ProgramParameters> tryLoadingProgramParameters(int argc, char** argv)
{
    auto optionalUserParameters = tryParsingUserParameters(argc, argv);
    if (!optionalUserParameters)
    {
        return boost::optional<ProgramParameters>();
    }

    const auto& userParams = *optionalUserParameters;
    assertValidity(userParams);

    InputPaths inputPaths(userParams.htsFilePath, userParams.referencePath, userParams.catalogPath);
    const string vcfPath = userParams.outputPrefix + ".vcf";
    const string jsonPath = userParams.outputPrefix + ".json";
    const string bamletPath = userParams.outputPrefix + "_realigned.bam";
    OutputPaths outputPaths(vcfPath, jsonPath, bamletPath);
    SampleParameters sampleParameters = decodeSampleParameters(userParams);
    HeuristicParameters heuristicParameters(
        userParams.regionExtensionLength, userParams.qualityCutoffForGoodBaseCall, userParams.skipUnaligned,
        userParams.alignerType);

    LogLevel logLevel;
    try
    {
        logLevel = decodeLogLevel(userParams.logLevel);
    }
    catch (std::logic_error)
    {
        const string message = "Log level must be set to either trace, debug, info, warn, or error";
        throw std::invalid_argument(message);
    }

    AnalysisMode analysisMode;
    try
    {
        analysisMode = decodeAnalysisMode(userParams.analysisMode);
    }
    catch (std::logic_error)
    {
        const string message = "Analysis mode must be set to either streaming or seeking";
        throw std::invalid_argument(message);
    }

    return ProgramParameters(inputPaths, outputPaths, sampleParameters, heuristicParameters, analysisMode, logLevel);
} */

/*
static void assertWritablePath(const string& pathEncoding)
{
    const fs::path path(pathEncoding);
    const fs::path pathToDirectory = path.parent_path();

    const bool thereIsNoDirectory = pathToDirectory.empty();
    const bool pathLeadsToExistingDirectory = fs::is_directory(pathToDirectory);
    const bool filenameIsValid = fs::portable_posix_name(path.filename().string());

    if (!filenameIsValid || (!thereIsNoDirectory && !pathLeadsToExistingDirectory))
    {
        throw std::invalid_argument(pathEncoding + " is not a valid output path");
    }
}*/

static void assertPathToExistingFile(const string& pathEncoding)
{
    const fs::path path(pathEncoding);
    const bool isPathToExistingFile = fs::exists(path) && !fs::is_directory(path);

    if (!isPathToExistingFile)
    {
        throw std::invalid_argument(pathEncoding + " is not a path to an existing file");
    }
}

void assertValidity(const MergeWorkflowParameters& parameters)
{
    assertPathToExistingFile(parameters.pathToReference());
    assertPathToExistingFile(parameters.pathToManifest());
}
