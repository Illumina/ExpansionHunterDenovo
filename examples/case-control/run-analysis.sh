#!/bin/bash


# Generate STR profiles for each sample
for bamlet in bamlets/*.bam
do
  sample=$(basename $bamlet)
  sample=${sample%.bam}

  ../../build/ExpansionHunterDenovo profile \
    --reads $bamlet \
    --reference reference.fasta \
    --output-prefix str-profiles/${sample}
done

# Merge STR profiles into multi-sample STR profile
../../build/ExpansionHunterDenovo merge \
  --reference reference.fasta \
  --manifest manifest.tsv \
  --output-prefix example_dataset


# Perform locus-based case-control comparison
../../scripts/casecontrol.py locus \
    --manifest manifest.tsv \
    --multisample-profile example_dataset.multisample_profile.json \
    --output example_dataset.casecontrol_locus.tsv


# Perform motif-based case-control comparison
../../scripts/casecontrol.py motif \
    --manifest manifest.tsv \
    --multisample-profile example_dataset.multisample_profile.json \
    --output example_dataset.casecontrol_motif.tsv

