# Merging single-sample profiles into multisample profiles

Once STR profiles have been computed for each sample in the dataset, they must
be merged together into a **multisample STR profile**. Multisample STR profiles
are produced by `ExpansionHunterDenovo merge` command which has the following
parameters.

| Required parameter | Description                                    |
|--------------------|------------------------------------------------|
| --reference        | The FASTA file to which the reads were aligned |
| --manifest         | TSV file with describing each sample           |
| --output-prefix    | Common prefix for the output files             |

| Optional parameter | Description                      | Default |
|--------------------|----------------------------------|:-------:|
| --min-unit-len     | Shortest repeat unit to consider | 2       |
| --max-unit-len     | Longest repeat unit to consider  | 20      |

## Manifest files

The manifest file is a tab-delimited file whose columns contain
sample id, case/control status, and absolute path to the EH Denovo
output file. For example,

```
sample1	case	/path/to/sample1.str_profile.json
sample2	case	/path/to/sample2.str_profile.json
sample3	case	/path/to/sample3.str_profile.json
sample4	control	/path/to/sample4.str_profile.json
sample5	control	/path/to/sample5.str_profile.json
sample6	control	/path/to/sample6.str_profile.json
```

## Multisample STR profiles

A multisample STR profile is a result of merging multiple single-sample STR
profiles. It is a JSON object containing counts of in-repeat reads across
all samples for each repeat motif. When two samples have overlapping or
nearby regions with anchored in-repeat reads corresponding to the same repeat
motif, these regions are merged together. By default, `merge` operation merges
regions located within 500bp of one another.

## Example

Let's assume we have two (single-sample) STR profiles that we want to combine.
The first profile is the file `sample1.str_profile.json` containing
```json
{
    "CCG": {
        "AnchoredIrrCount": 7,
        "IrrPairCount": 3,
        "RegionsWithIrrAnchors": {
            "chr1:100-1000": 7
        },
        "RepeatUnit": "CCG",
        "ReadLength": 150,
        "Depth": 30.0
    }
}
```

And second profile is the file `sample2.str_profile.json` containing

```json
{
    "AGG": {
        "AnchoredIrrCount": 7,
        "IrrPairCount": 2,
        "RegionsWithIrrAnchors": {
            "chr2:5000-5300": 7
        },
        "RepeatUnit": "AGG"
    },
    "CCG": {
        "AnchoredIrrCount": 3,
        "IrrPairCount": 10,
        "RegionsWithIrrAnchors": {
            "chr1:1150-1200": 3
        },
        "RepeatUnit": "CCG",
        "ReadLength": 150,
        "Depth": 40.0
    }
}
```

To merge these profiles together we create a manifest file `manifest.txt`

```
sample1	case	/path/to/sample1.str_profile.json
sample2	control	/path/to/sample2.str_profile.json
```

and then run

```bash
ExpansionHunterDenovo merge \
  --reference reference.fa \
  --manifest manifest.txt \
  --output-prefix example_dataset
```

The program generates a file `example_dataset.multisample_str_profile.json`
with the following content

```json
{
    "Counts": {
        "AGG": {
            "IrrPairCounts": {
                "sample2": 2
            },
            "RegionsWithIrrAnchors": {
                "chr2:5000-5300": {
                    "sample2": 7
                }
            }
        },
        "CCG": {
            "IrrPairCounts": {
                "sample1": 3,
                "sample2": 10
            },
            "RegionsWithIrrAnchors": {
                "chr1:100-1200": {
                    "sample1": 7,
                    "sample2": 3
                },
            }
        }
    },
    "Parameters": {
        "Depths": {
            "sample1": 30.0,
            "sample2": 40.0
        },
        "ReadLengths": {
            "sample1": 150,
            "sample2": 150
        }
    }
}
```