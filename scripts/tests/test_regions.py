#
# ExpansionHunter Denovo
# Copyright 2016-2019 Illumina, Inc.
# All rights reserved.
#
# Author: Egor Dolzhenko <edolzhenko@illumina.com>,
#         Michael Eberle <meberle@illumina.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#

import sys
import pytest

from core.regiontools import *


class TestFeatureCounts(object):
    def test_basic_operations(self):
        counts = FeatureCounts({"S1": 5, "S2": 10})
        assert len(counts) == 2
        assert counts["S1"] == 5 and counts["S2"] == 10

    def test_repr(self):
        counts = FeatureCounts({"S1": 5, "S2": 10})
        assert repr(counts) == "S1=5,S2=10"

    def test_dictionary_conversion(self):
        feature_counts = FeatureCounts({"S1": 5, "S2": 10})
        assert dict(feature_counts) == {"S1": 5, "S2": 10}

    def test_extension_with_disjoint_samples(self):
        counts_a = FeatureCounts({"S1": 5, "S2": 10})
        counts_b = FeatureCounts({"S3": 4, "S4": 10})
        counts_a.combine(counts_b)
        assert counts_a == FeatureCounts({"S1": 5, "S2": 10, "S3": 4, "S4": 10})

    def test_extension_with_overlapping_samples(self):
        counts_a = FeatureCounts({"S1": 5, "S2": 10})
        counts_b = FeatureCounts({"S2": 5, "S3": 20})
        counts_a.combine(counts_b)
        assert counts_a == FeatureCounts({"S1": 5, "S2": 15, "S3": 20})


class TestRegionProperties(object):
    def test_unset_by_default(self):
        region = Region()
        assert region.is_unset()

    def test_representation(self):
        region = Region("chr1", 10, 20)
        assert repr(region) == "chr1:10-20"

    def test_representation_with_feature_counts(self):
        feature_counts = FeatureCounts({"Sample1": 5, "Sample2": 10})
        region = Region("chr1", 10, 20, feature_counts)
        assert repr(region) == "chr1:10-20(Sample1=5,Sample2=10)"

    def test_comparison(self):
        assert Region("chr1", 10, 20) < Region("chr2", 10, 20)
        assert Region("chr1", 10, 20) < Region("chr1", 15, 20)
        assert Region("chr1", 10, 15) < Region("chr1", 10, 20)

    def test_equality(self):
        assert Region("chr1", 10, 20) == Region("chr1", 10, 20)
        assert Region("chr1", 10, 20, FeatureCounts({"S1": 1})) == Region(
            "chr1", 10, 20, FeatureCounts({"S1": 1})
        )
        assert Region("chr1", 10, 20, FeatureCounts({"S1": 1})) != Region(
            "chr1", 10, 20, FeatureCounts({"S1": 2})
        )

    def test_initializes_from_denovo_record(self):
        region = create_region_from_denovo_record("S1", ("10:100-104", 1))
        assert region == Region("10", 100, 104, FeatureCounts({"S1": 1}))


class TestDistanceCalculation(object):
    def test_overlapping_regions(self):
        region_a = Region("chr1", 1, 10)
        region_b = Region("chr1", 5, 15)
        assert compute_distance(region_a, region_b) == 0

    def test_disjoint_regions(self):
        region_a = Region("chr1", 50, 70)
        region_b = Region("chr1", 0, 20)
        assert compute_distance(region_a, region_b) == 30
        assert compute_distance(region_b, region_a) == 30

    def test_regions_on_different_chroms(self):
        region_a = Region("chr1", 10, 20)
        region_b = Region("chr2", 10, 20)
        assert compute_distance(region_a, region_b) == sys.maxsize

    def test_unset_region_rases_exception(self):
        region_a = Region("chr1", 10, 20)
        region_b = Region()
        with pytest.raises(Exception) as excinfo:
            compute_distance(region_a, region_b)
        msg = str(excinfo.value)
        assert msg == "Cannot compute distance between unset regions"


class TestRegionCollection(object):
    def test_basic_operations(self):
        regions = RegionCollection([Region("1", 1, 2), Region("2", 3, 4)])
        assert len(regions) == 2
        assert regions[0] == Region("1", 1, 2)
        assert regions[1] == Region("2", 3, 4)

    def test_has_no_elements_by_default(self):
        assert len(RegionCollection()) == 0

    def test_dictionary_conversion(self):
        regions = RegionCollection(
            [
                Region("1", 0, 10, FeatureCounts({"S1": 5, "S2": 1})),
                Region("1", 200, 250, FeatureCounts({"S1": 10})),
            ]
        )
        assert regions.as_dict() == {
            "1:0-10": {"S1": 5, "S2": 1},
            "1:200-250": {"S1": 10},
        }

    def test_extend(self):
        regions_a = RegionCollection(
            [
                Region("1", 0, 10, FeatureCounts({"S1": 5, "S2": 1})),
                Region("1", 200, 250, FeatureCounts({"S1": 10})),
            ]
        )

        regions_b = RegionCollection([Region("1", 0, 5, FeatureCounts({"S1": 7}))])
        regions_a.extend(regions_b)

        expected_regions = RegionCollection(
            [
                Region("1", 0, 10, FeatureCounts({"S1": 5, "S2": 1})),
                Region("1", 200, 250, FeatureCounts({"S1": 10})),
                Region("1", 0, 5, FeatureCounts({"S1": 7})),
            ]
        )

        assert regions_a == expected_regions

    def test_create_from_denovo_record(self):
        regions = create_region_collection_from_denovo_record(
            "S1", {"1:1-10": 10, "2:1-20": 20}
        )

        expected_regions = RegionCollection(
            [
                Region("1", 1, 10, FeatureCounts({"S1": 10})),
                Region("2", 1, 20, FeatureCounts({"S1": 20})),
            ]
        )
        assert regions == expected_regions


class TestMerging(object):
    def test_sorted_overlapping_regions(self):
        regions = RegionCollection(
            [
                Region("chr1", 10, 20, FeatureCounts({"S1": 1})),
                Region("chr1", 15, 25, FeatureCounts({"S2": 2})),
                Region("chr1", 20, 35, FeatureCounts({"S3": 3})),
            ]
        )
        regions.merge()

        expected_regions = RegionCollection(
            [Region("chr1", 10, 35, FeatureCounts({"S1": 1, "S2": 2, "S3": 3}))]
        )
        assert regions == expected_regions

    def test_unsorted_overlapping_regions(self):
        regions = RegionCollection(
            [
                Region("chr1", 15, 25, FeatureCounts({"S1": 1})),
                Region("chr1", 10, 20, FeatureCounts({"S1": 1})),
                Region("chr1", 20, 35, FeatureCounts({"S1": 1})),
            ]
        )
        regions.merge()

        expected_regions = RegionCollection(
            [Region("chr1", 10, 35, FeatureCounts({"S1": 3}))]
        )
        assert regions == expected_regions

    def test_disjoint_regions(self):
        regions = RegionCollection(
            [
                Region("chr1", 15, 25, FeatureCounts({"S1": 1, "S2": 1})),
                Region("chr2", 10, 20, FeatureCounts({"S1": 1})),
                Region("chr1", 200, 350, FeatureCounts({"S1": 1, "S2": 1})),
            ]
        )
        regions.merge()

        expected_regions = RegionCollection(
            [
                Region("chr1", 15, 350, FeatureCounts({"S1": 2, "S2": 2})),
                Region("chr2", 10, 20, FeatureCounts({"S1": 1})),
            ]
        )
        assert regions == expected_regions

    def test_multiple_regions(self):
        regions = RegionCollection(
            [
                Region("1", 200, 250, FeatureCounts({"S1": 10})),
                Region("1", 500, 550, FeatureCounts({"S1": 20})),
                Region("1", 0, 10, FeatureCounts({"S1": 5, "S2": 1})),
                Region("1", 1100, 1200, FeatureCounts({"S2": 30})),
                Region("2", 1100, 1200, FeatureCounts({"S2": 20})),
            ]
        )
        regions.merge(max_dist=300)
        expected_regions = RegionCollection(
            [
                Region("1", 0, 550, FeatureCounts({"S1": 35, "S2": 1})),
                Region("1", 1100, 1200, FeatureCounts({"S2": 30})),
                Region("2", 1100, 1200, FeatureCounts({"S2": 20})),
            ]
        )
        assert regions == expected_regions

    def test_included_regions(self):
        regions = RegionCollection(
            [
                Region("1", 100, 200, FeatureCounts({"S1": 1})),
                Region("1", 90, 300, FeatureCounts({"S1": 1})),
            ]
        )
        regions.merge()
        expected_regions = RegionCollection(
            [Region("1", 90, 300, FeatureCounts({"S1": 2}))]
        )
        assert regions == expected_regions

