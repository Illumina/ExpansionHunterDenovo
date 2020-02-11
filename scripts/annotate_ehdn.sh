#!/usr/bin/env bash
#
# ExpansionHunter Denovo results annotation
#
# Author: Mark Bennett <mark.bennett@wehi.edu.au>
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

usage() {
  echo -n "annotate_ehdn.sh [options]
  
Available options:
  -h | --help                   Print help message
  --ehdn-results                ExpansionHunterDenovo secondary locus analysis tsv
  --ehdn-annotated-results      ExpansionHunterDenovo annotated output filename
  --annovar-annotate-variation  ANNOVAR annotate_variation.pl script path
  --annovar-humandb             ANNOVAR humandb directory
  --annovar-buildver            ANNOVAR buildver option (hg19, hg38, ...)

Annotation script requires ANNOVAR installation.
1. Register and download ANNOVAR from http://annovar.openbioinformatics.org
2. Download refGene annotation database to 'humandb/' such as:
   'annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/'

"
}

# Process command line options

PROVIDE_EHDN_RESULTS=false
PROVIDE_EHDN_ANNOTATED=false
PROVIDE_ANNOVAR_ANNOTATE_VARIATION=false
PROVIDE_ANNOVAR_HUMANDB=false
PROVIDE_ANNOVAR_BUILDVER=false

while [[ $1 = -?* ]]
do
  case $1 in
    -h|--help)
        usage >&2
        exit 0 ;;
    --ehdn-results) 
        shift
        EHDN_RESULTS=${1}
        PROVIDE_EHDN_RESULTS=true ;;
    --ehdn-annotated-results)
        shift
        EHDN_ANNOTATED=${1}
        PROVIDE_EHDN_ANNOTATED=true ;;
    --annovar-annotate-variation)
        shift
        ANNOVAR_ANNOTATE_VARIATION=${1}
        PROVIDE_ANNOVAR_ANNOTATE_VARIATION=true ;;
    --annovar-humandb)
        shift
        ANNOVAR_HUMANDB=${1}
        PROVIDE_ANNOVAR_HUMANDB=true ;;
    --annovar-buildver)
        shift
        ANNOVAR_BUILDVER=${1}
        PROVIDE_ANNOVAR_BUILDVER=true ;;
    *) 
        echo -e "Error: invalid option '$1'\nTry 'annotate_ehdn.sh --help' for information\n"  >&2
        exit 1 ;;
  esac
  shift
done

# Check all required argements specified
if [[ ("$PROVIDE_EHDN_RESULTS" == true) && ("$PROVIDE_EHDN_ANNOTATED" == true) && ("$PROVIDE_ANNOVAR_ANNOTATE_VARIATION" == true) && ("$PROVIDE_ANNOVAR_HUMANDB" == true) && ("$PROVIDE_ANNOVAR_BUILDVER" == true) ]]
then
    :
else
    echo -e "Error: missing required argument(s)\nTry 'annotate_ehdn.sh --help' for information\n"  >&2
    exit 1
fi

# Check that EHdn results header is compatible (eg 'motif' analysis invalid)
# (case-control and outlier output have inconsistent headers)
HEADER_COL2=$(head -n1 $EHDN_RESULTS | cut -f 2)
HEADER_COL3=$(head -n1 $EHDN_RESULTS | cut -f 3)
if [ "$HEADER_COL2" != "start" ] || [ "$HEADER_COL3" != "end" ]
then
    echo -e "Error: unexpected column names. EHdn result file does not appear to have the correct format\n" >&2
    exit 1
fi
    

AVINPUT=${EHDN_ANNOTATED}.temp${RANDOM}.avinput

# Convert to ANNOVAR input format: add 'ref' and 'alt' columns containing '0'
tail -n +2 $EHDN_RESULTS | awk '{OFS=FS="\t"} $4="0\t0\t"$4' > $AVINPUT

# Run ANNOVAR
$ANNOVAR_ANNOTATE_VARIATION \
    -geneanno -dbtype refGene -buildver $ANNOVAR_BUILDVER \
    $AVINPUT $ANNOVAR_HUMANDB > /dev/null 2>&1

# Reformat output: move gene and region back and remove temp ref/alt columns
head -n 1 $EHDN_RESULTS | awk '{OFS=FS="\t"} $5="gene\tregion\t"$5' \
    > $EHDN_ANNOTATED
awk 'BEGIN {OFS="\t"}; {print $3, $4, $5, $8, $2, $1, $9, $10, $11}' \
    $AVINPUT.variant_function >> $EHDN_ANNOTATED

# Clean up
rm $AVINPUT
rm $AVINPUT.exonic_variant_function
rm $AVINPUT.log
rm $AVINPUT.variant_function


