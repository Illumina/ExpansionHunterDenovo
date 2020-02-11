# Before you begin

ExpansionHunter Denovo (EHdn) can be used to analyze a collection of BAM/CRAM
files containing alignments of short (100-200bp) reads. For best results, the
samples should be sequenced on the same instrument to similar coverage of at
least 30x. All data should be aligned with the **same short-read aligner**
ideally without any post-processing steps such as indel realignment or
recalibration.

Your dataset should contain one or more samples that are suspected to harbor a
repeat expansion and a set of controls. If the controls are not available,
consider using
[Illumina Polaris](https://github.com/Illumina/Polaris/wiki/HiSeqX-Diversity-Cohort)
dataset.
