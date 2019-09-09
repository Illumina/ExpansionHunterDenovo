# Case-control analysis



The case-control analysis is performed by a Python3 script
`case_control_analysis.py` located inside `scripts` directory. The script can
be run like so:

```bash
python3 case_control_analysis.py \
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
| --test-method      | Method of calculating Wilcoxon Rank-Sum Test p-value         | normal  |
| --num-resamples    | Number of iterations for the resampling test                 | 1000000 |

## Outputs

The program produces two output files. One of them summarizes per-locus
comparison of in-repeat reads. The other file summarizes the overall genome-wide
comparison of motifs.
