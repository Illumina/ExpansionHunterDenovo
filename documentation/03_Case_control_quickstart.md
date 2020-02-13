# Quickstart: Case-control analysis

A case-control analysis is appropriate when a significant subset of cases is
expected to contain expansions of the same repeat. For example, if a
case-control analysis is applied to a dataset consisting of ALS patients and
healthy controls, it is expected to flag the GGCCCC repeat in *C9orf72* gene
as significant.

Let's perform a case-control analysis of a tiny simulated dataset consisting
of seven samples. The data we will be using is contained inside directory
`ExpansionHunterDenovo/examples/case-control`. We assume that this is our
working directory.

We start by generating an **STR profile** for the first sample:

```bash
/path/to/ExpansionHunterDenovo profile \
        --reads bamlets/sample1.bam \
        --reference reference.fasta \
        --output-prefix str-profiles/sample1 \
        --min-anchor-mapq 50 \
        --max-irr-mapq 40
```

This command will produce an STR profile `sample1.str_profile.json` inside the
`str-profiles` directory containing information about STR longer than the read
length in `sample1`. STR profiles for the remaining six BAM files can be
computed the same way.

We next merge the STR profiles together into a **multi-sample STR profile**.
This requires us to create the manifest file `manifest.tsv` describing the
dataset. The manifest file contains columns for sample identifier,
case-control status, and path to the associated STR profile for each sample.

```bash
sample1 case    str-profiles/sample1.str_profile.json
sample2 case    str-profiles/sample2.str_profile.json
sample3 case    str-profiles/sample3.str_profile.json
sample4 control str-profiles/sample4.str_profile.json
sample5 control str-profiles/sample5.str_profile.json
sample6 control str-profiles/sample6.str_profile.json
sample7 control str-profiles/sample7.str_profile.json
```

The merge operation is performed like so:

```bash
/path/to/ExpansionHunterDenovo merge \
        --reference reference.fasta \
        --manifest manifest.tsv \
        --output-prefix example_dataset
```

The resulting multisample STR profile is named
`example_dataset.multisample_profile.json`. It contains information
about long repeats across the entire dataset and sample-specific
parameters like read length and read coverage depth.

The next step is to actually compare STR lengths between cases and controls. Two
types of comparisons are supported: locus-based comparison and motif-based
comparison. These comparisons are performed by a Python3 script `casecontrol.py`
located in the `scripts/` directory.

- The locus-based comparison can distinguish between repeats longer
than the read length but shorter than the fragment length. It can be
run like so:

    ```bash
    /path/to/scripts/casecontrol.py locus \
            --manifest manifest.tsv \
            --multisample-profile example_dataset.multisample_profile.json \
            --output example_dataset.casecontrol_locus.tsv \
            [--target-regions target_regions.bed]
    ```

  where `--target-regions` is an optional parameter for restricting the analysis
  to regions in the provided a BED file.

  The resulting file `example_dataset.casecontrol_locus.tsv` contains the
  following information.
  
  | Column      | Description                         |
  |-------------|-------------------------------------|
  | contig      | Contig of the repeat region         |
  | start       | Approximate start of the repeat     |
  | end         | Approximate end of the repeat       |
  | motif       | Inferred repeat motif               |
  | pvalue      | P-value from Wilcoxon rank-sum test |
  | bonf_pvalue | P-value after Bonferroni correction |
  | counts      | Depth-normalized counts of anchored in-repeat reads for each sample (omitting samples with zero count) |

- The motif-based comparison analyzes the overall enrichment of
genomes with repeats longer than the fragment length. It can be run like so:

    ```bash
    path/to/scripts/casecontrol.py motif \
            --manifest manifest.tsv \
            --multisample-profile example_dataset.multisample_profile.json \
            --output example_dataset.casecontrol_motif.tsv
    ```

  The resulting file `example_dataset.casecontrol_motif.tsv` contains the
  following information.
  
  | Column      | Description                         |
  |-------------|-------------------------------------|
  | motif       | Inferred repeat motif               |
  | pvalue      | P-value from Wilcoxon rank-sum test |
  | bonf_pvalue | P-value after Bonferroni correction |
  | counts      | Depth-normalized counts of in-repeat read pairs for each sample (omitting samples with zero count) |
