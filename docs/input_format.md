# Input Format

This page describes the format of the GWAS summary statistics data (specified
using the `--gcor` flag) andthe reference panel required by X-LDSC (specified
using the `--bfile` format).

## GWAS summary statistics data

GWAS summary statistics data should be in plain text or gzipped text format
containing the following columns:

* SNP - rs ID of the SNP (e.g. rs62442).
* CHR - Chromosome number of the SNP. This should be a number between 1 and 22.
* BP - Base pair position of the SNP.
* A1 - Effect allele of the SNP. The sign of the Z-score is with respect to this allele.
* A2 - The other allele of the SNP.
* Z - The Z-score of the SNP.
* N - Sample size of the SNP.

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>:  This file format is compatible with <a href="https://github.com/bulik/ldsc">LDSC</a>.
Other columns [e.g. MAF, INFO, etc.] may be included in the file, but will not
be used. It is also recommended [although not required] that the summary data
files are sorted by their chromosomes and base pairs. All SNPs with either
duplicate ID or position will be removed before any analysis. )
</div>

## Reference panel

Reference panels should be in [PLINK format](https://www.cog-genomics.org/plink/2.0/input#bed)
(specified using the `--bfile` flag).

The following is a list of publicly available reference panels.

* [1000 Genomes Project](http://www.internationalgenome.org/data/)
* [UK10K](https://www.uk10k.org/data_access.html)

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: X-LDSC requires the centimorgan information of each SNP to
accurate estimation of LD scores. )
</div>
