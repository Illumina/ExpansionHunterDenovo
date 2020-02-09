#!/bin/bash

for bamlet in bamlets/*.bam
do
  sample=$(basename $bamlet)
  sample=${sample%.bam}

  ../../build/ExpansionHunterDenovo profile \
    --reads $bamlet \
    --reference reference.fasta \
    --output-prefix str-profiles/${sample}
done

../../build/ExpansionHunterDenovo merge \
  --reference reference.fasta \
  --manifest manifest.tsv \
  --output-prefix example_dataset

../../scripts/outlier.py locus \
  --manifest manifest.tsv \
  --multisample-profile example_dataset.multisample_profile.json \
  --output example_dataset.outlier_locus.tsv

../../scripts/outlier.py motif \
  --manifest manifest.tsv \
  --multisample-profile example_dataset.multisample_profile.json \
  --output example_dataset.outlier_motif.tsv

