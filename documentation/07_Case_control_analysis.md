# Case-control analysis

The case-control analysis is performed by a Python3 script `casecontrol.py`
located inside `scripts` directory. The locus-based analysis can be run
like so:

```bash
python3 case_control_analysis.py locus \
    --manifest manifest.txt \
    --multisample-profile multisample_profile.json \
    --output-prefix output
```

The command to run the motif-based analysis is nearly identical:

```bash
python3 case_control_analysis.py motif \
    --manifest manifest.txt \
    --multisample-profile multisample_profile.json \
    --output-prefix output
```

The input parameters manifest.txt and multisample_profile.json are as
[described previously](04_Merging_profiles.md). 

| Optional parameter | Description                                                  | Default |
|--------------------|--------------------------------------------------------------|:-------:|
| --min-count        | Minimum number reads in a region for downstream analysis     | 5       |
| --target-regions   | BED file with regions to which analysis should be restricted | NA      |
| --test-method      | Method of calculating Wilcoxon Rank-Sum Test p-value*         | normal  |

\* The default value `normal` invokes the Normal approximation appropriate for
larger samples. To compute the p-value directly for smaller samples, use
`permute_<N>` where N is the number of permutations. For example,
permute_1000000 invokes a test with 1000000 permutations.

## Outputs

The program produces two output files. One of them summarizes per-locus
comparison of in-repeat reads. The other file summarizes the overall genome-wide
comparison of motifs.
