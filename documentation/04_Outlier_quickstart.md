# Quickstart: Outlier analysis

If cases consist of samples from patients with diverse phenotypes, it might be
appropriate to assume that there is no enrichment for any specific expansion
and hence the case-control analysis is not appropriate. In this situation, an
outlier analysis can be used to flag repeats that are expanded in a small
proportion of cases compared to the rest of the dataset.

The example dataset used for this section is located in
`ExpansionHunterDenovo/examples/outlier`. We assume that individual STR
profiles and the multisample STR profile are computed as described in the
[quick start guide to case-control analysis](03_Case_control_quickstart.md).

The outlier analysis is performed by a Python3 script `outlier.py` located
in the `scripts/` directory. Two types of comparisons are supported:
locus-based comparison and motif-based comparison.

- The locus-based comparison can distinguish between repeats longer
than the read length but shorter than the fragment length. It can be
run like so:

    ```bash
    /path/to/scripts/outlier.py locus \
            --manifest manifest.tsv \
            --multisample-profile example_dataset.multisample_profile.json \
            --output example_dataset.outlier_locus.tsv
    ```

  The resulting file `example_dataset.outlier_locus.tsv` contains the
  following information.

  | Column           | Description                                                      |
  |------------------|------------------------------------------------------------------|
  | contig           | Contig of the repeat region                                      |
  | start            | Approximate start of the repeat                                  |
  | end              | Approximate end of the repeat                                    |
  | motif            | Inferred repeat motif                                            |
  | top_case_zscore  | Top z-score of a case sample                                     |
  | high_case_counts | Counts of case samples corresponding to z-score greater than 1.0 |
  | counts           | Nonzero counts for all samples                                   |

- The motif-based comparison analyzes the overall enrichment of genomes with
repeats longer than the fragment length. It can be run like so:

    ```bash
    path/to/scripts/outlier.py motif \
            --manifest manifest.tsv \
            --multisample-profile example_dataset.multisample_profile.json \
            --output example_dataset.outlier_motif.tsv \
            [--target-regions target_regions.bed]
    ```

  where `--target-regions` is an optional parameter for restricting the analysis
  to regions in the provided a BED file.

  The resulting file `example_dataset.outlier_motif.tsv` contains the following
  information.

  | Column           | Description                                                      |
  |------------------|------------------------------------------------------------------|
  | motif            | Inferred repeat motif                                            |
  | top_case_zscore  | Top z-score of a case sample                                     |
  | high_case_counts | Counts of case samples corresponding to z-score greater than 1.0 |
  | counts           | Nonzero counts for all samples                                   |
