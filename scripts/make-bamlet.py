#!/usr/bin/env python3
#
# ExpansionHunter Denovo
# Copyright 2016-2020 Illumina, Inc.
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

import os
import sys
import pysam
import argparse

def is_close(chrom, pos, region):
    max_dist = 1000
    reg_chrom, reg_start, reg_end = region

    if chrom != reg_chrom:
        return False

    dist = min(abs(pos - reg_start), abs(pos - reg_end))
    if dist > max_dist:
        return False

    return True

def jump_for_mate(bam_path, chrom, pos, read_name):
    bam = pysam.AlignmentFile(bam_path, 'rb')   
    for al in bam.fetch(chrom, pos, pos + 1):
        if al.is_secondary:
            continue
        if al.query_name == read_name:
            bam.close()
            return al
    bam.close()
    print('[WARNING: Could not locate {} at {}:{}]'.format(read_name, chrom, pos))
    return False 

def extract_region(region, bam_path, bamlet_path):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    bamlet = pysam.AlignmentFile(bamlet_path, 'wb', template=bam)

    mates = {}

    for al in bam.fetch(*region):
        if al.is_secondary:
             continue

        read_name = al.query_name
        if read_name not in mates:
            mates[read_name] = []
        mates[read_name].append(al)
        #assert len(mates[read_name]) <= 2, 'ERROR: Found more than two reads named ' + read_name

    for read_name in mates:
        reads = mates[read_name]
        if len(reads) == 2:
            continue
        read = reads[0] 
        mate_chrom = read.next_reference_name
        mate_pos = read.next_reference_start

        if not is_close(mate_chrom, mate_pos, region):
            print('[Looking for mate of {} in {}:{}]'.format(read_name, mate_chrom, mate_pos))
            mate = jump_for_mate(bam_path, mate_chrom, mate_pos, read_name)
            if mate:
                reads.append(mate)

    for read_name, rec in mates.items():
        #assert len(rec) <= 2, 'ERROR: Found more than two reads named ' + read_name
        for read in rec:
            bamlet.write(read)

    bam.close()
    bamlet.close()

parser = argparse.ArgumentParser(description='A script to generate BAMlets')
parser.add_argument('--bam', type=str, nargs=1,
                    required=True, help='Input BAM file')
parser.add_argument('--region', type=str, nargs=1,
                    required=True, help='Region from which to extract reads (chr:start-end)')
parser.add_argument('--bamlet', type=str, nargs=1,
                    required=True, help='Output BAMlet')

args = parser.parse_args()
bam_path = args.bam[0]
region = args.region[0]
bamlet_path = args.bamlet[0]

chrom, start, end = region.replace(':', ' ').replace('-', ' ').split()
extension_len = 2000
start, end = int(start) - extension_len, int(end) + extension_len
region = (chrom, start, end)

extract_region(region, bam_path, bamlet_path)
