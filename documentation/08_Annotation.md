# Annotate locus-based analysis results

This guide provide an example of one method by which the ExpansionHunter Denovo
case-control or outlier locus-based analysis results can be annotated with
information about the gene and genomic region for each repeat locus identified.

A script to annotate ExpansionHunter Denovo locus-based results using
[ANNOVAR](http://annovar.openbioinformatics.org) is provided here
`ExpansionHunterDenovo/scripts/annotate_ehdn.sh`.

**Note:** Genomic regions reported by this annotation script need to be
interpreted with caution.

Regions identified by ExpansionHunter Denovo are defined based on the
locations of the supporting anchored in-repeat reads and do not exactly
correspond to location of the repeat expansion itself. ANNOVAR denotes
a region as 'exonic' if any part of the region overlaps with an exon.

A locus annotated as 'exonic' does not guarantee that the repeat expansion
itself overlaps a coding region, only that the region defined by the
supporting anchored in-repeats reads overlaps with an exon. However, it does
suggest the expansion lies very close to a coding region.

## Setup ANNOVAR

The provided annotation script requires an existing ANNOVAR installation and the
refGene annotation database to be downloaded.

1. [Register and download ANNOVAR](http://download.openbioinformatics.org/annovar_download_form.php)
2. Unpack ANNOVAR in your home directory or any other directory in your `PATH`

    ```bash
    tar xvfz annovar.latest.tar.gz
    mv annovar ~/annovar
    ```

3. Download the refGene annotation database for the genome build you wish to use:

    ```bash
    perl ~/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ~/annovar/humandb/
    perl ~/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene ~/annovar/humandb/
    ```

4. The remainder of this guide assumes that ANNOVAR is installed in the
   directory `~/annovar/` and the annotation databases are downloaded to
   `~/annovar/humandb/`

## Usage

ExpansionHunter Denovo [case-control](03_Case_control_quickstart.md) or
[outlier](04_Outlier_quickstart.md) locus-based analysis is performed as
described elsewhere in [documentation](00_Introduction.md).

The following command annotates an EHdn results file called
`locus_output.tsv` to create an annotated output file
`locus_output_annotated.tsv`.

```bash
bash annotate_ehdn.sh \
    --ehdn-results locus_output.tsv \
    --ehdn-annotated-results locus_output_annotated.tsv \
    --annovar-annotate-variation ~/annovar/annotate_variation.pl \
    --annovar-humandb ~/annovar/humandb \
    --annovar-buildver hg19
```

Here is a description of the required arguments.

| Parameter                  | Description                                                          |
|----------------------------|----------------------------------------------------------------------|
| ehdn-results               | Input filename of EHdn locus analysis results                        |
| ehdn-annotated-results     | Output filename for annotated EHdn locus analysis results            |
| annovar-annotate-variation | location of annotate_variation.pl ANNOVAR script                     |
| annovar-humandb            | location of directory containing ANNOVAR refGene annotation database |
| annovar-builder            | reference genome build used by ANNOVAR ('hg18', 'hg19' or 'hg38')    |
