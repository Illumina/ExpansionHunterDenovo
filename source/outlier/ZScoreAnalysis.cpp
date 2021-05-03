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

#include "outlier/ZScoreAnalysis.hh"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

using std::string;
using std::vector;

ZScoreAnalysisResults analyzeZScores(const Manifest& manifest, const AnchoredIrrCounts& irrCounts)
{
    vector<double> rawCounts;
    rawCounts.reserve(manifest.size());

    for (const auto& sampleAndStatus : manifest)
    {
        const auto& sample = sampleAndStatus.first;
        auto record = irrCounts.countBySample.find(sample);
        if (record != irrCounts.countBySample.end())
        {
            rawCounts.push_back(record->second);
        }
        else
        {
            rawCounts.push_back(0);
        }
    }

    const int kNumResamples = 100;
    vector<double> resample;
    resample.reserve(rawCounts.size());
    const int targetQuantileIndex = static_cast<int>(std::round((rawCounts.size() - 1) * 0.95));

    std::random_device device;
    std::mt19937 engine { device() };
    std::uniform_int_distribution<int> randomIndex(0, rawCounts.size() - 1);

    vector<double> quantiles;
    quantiles.reserve(kNumResamples);
    for (int resampleIndex = 0; resampleIndex != kNumResamples; ++resampleIndex)
    {
        resample.clear();
        for (int countIndex = 0; countIndex != rawCounts.size(); ++countIndex)
        {
            resample.push_back(rawCounts[randomIndex(engine)]);
        }

        std::nth_element(resample.begin(), resample.begin() + targetQuantileIndex, resample.end());
        quantiles.push_back(resample[targetQuantileIndex]);
    }

    // std::sort(resample.begin(), resample.end());
    for (auto value : quantiles)
    {
        std::cerr << value << " ";
    }
    std::cerr << std::endl;

    const double mean = std::accumulate(quantiles.begin(), quantiles.end(), 0.0) / quantiles.size();

    double variance = 0;
    for (auto quantile : quantiles)
    {
        variance += std::pow(quantile - mean, 2);
    }
    variance /= quantiles.size();
    const double sigma = std::max(std::sqrt(variance), 1.0);

    std::cerr << mean << " " << sigma << std::endl;

    double topZScore = -1;
    vector<string> casesWithHighCounts;
    for (const auto& sampleAndCount : irrCounts.countBySample)
    {
        const auto& sample = sampleAndCount.first;
        const double count = sampleAndCount.second;
        const auto sampleStatus = manifest.find(sample)->second.status;
        if (sampleStatus == SampleStatus::kCase)
        {
            const double zScore = (count - mean) / sigma;
            if (zScore > 1.0)
            {
                topZScore = std::max(topZScore, zScore);
                casesWithHighCounts.push_back(sample);
            }
        }
    }

    return { topZScore, casesWithHighCounts };
}
