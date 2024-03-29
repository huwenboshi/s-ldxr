# Input Format

This page describes the format of GWAS summary statistics data (specified
by the `--gcor` flag), reference panel (specified by the `--bfile` flag),
annotation files (specified by the `--annot` flag), LD score files (
specified by the `--ref-ld-chr` and `--w-ld-chr` flag), and minor allele
frequency files (specified by the `--frqfile` flag).

## GWAS summary statistics data

GWAS summary statistics data should be in plain text or gzipped text format
containing the following columns:

* `SNP` rs ID of the SNP (e.g. rs62442).
* `CHR` Chromosome number of the SNP. This should be a number between 1 and 22.
* `BP` Base pair position of the SNP.
* `A1` Effect allele of the SNP. The sign of the Z-score is with respect to this allele.
* `A2` The other allele of the SNP.
* `Z` The Z-score of the SNP.
* `N` Sample size of the SNP.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>:  This file format is compatible with <a href="https://github.com/bulik/ldsc">LDSC</a>.
Other columns (e.g. MAF, INFO, etc.) may be included in the file, but will not
be used. It is also recommended [although not required] that the summary data
files are sorted by their chromosomes and base pairs. All SNPs with either
duplicate ID or position will be removed before any analysis.
</div>

<div style="background-color:rgba(240, 128, 128, 0.2);">
( <b>Note</b>: It is important that the alleles of each SNP in the GWAS
summary statistics file are consistent across the two populations.)
</div>

## Reference panel

Reference panels should be in [PLINK format](https://www.cog-genomics.org/plink/2.0/input#bed)
(specified using the `--bfile` flag).

The following is a list of publicly available reference panels.

* [1000 Genomes Project](http://www.internationalgenome.org/data/)
* [UK10K](https://www.uk10k.org/data_access.html)

1000 Genomes Project reference genotype data for East Asian and European
populations can be downloaded [here](https://data.broadinstitute.org/alkesgroup/S-LDXR/1000G_EAS_EUR.tar.gz).

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: S-LDXR requires the centimorgan information of each SNP to
accurate estimation of LD scores.
</div>

## Annotation files

Annotation files should be in gzipped text format. The first 4 columns should
always be:

* `CHR` Chromosome number of the SNP.
* `BP` Base pair position of the SNP.
* `SNP` rs ID of the SNP (e.g. rs19800731).
* `CM` Position of the SNP in centimogan.

The first annotation (usually the base annotation) should start on the 5th
column.

baseline-LD-X model annotations can be downloaded
[here](https://data.broadinstitute.org/alkesgroup/S-LDXR/baseline-LD-X.tar.gz).

## LD score files

LD score files should be in gzipped text format. The first 3 columns should
always be:

* `CHR` Chromosome number of the SNP.
* `SNP` rs ID of the SNP (e.g. rs19790919).
* `BP` Base pair position of the SNP.

LD scores of the first annotation (usually the base annotation) should start
on the 4th column.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: The regression weight file (specified by the --w-ld-chr flag)
has the same format as LD score files.
</div>

LD scores for the baseline-LD-X model annotations can be downloaded
[here](https://data.broadinstitute.org/alkesgroup/S-LDXR/baseline-LD-X.tar.gz).

## Minor allele frequency files

Minor allele frequency files format is the same as that of PLINK. The columns
are:

* `CHR` Chromosome number of the SNP.
* `SNP` rs ID of the SNP (e.g. rs19800301).
* `A1` Non-effect allele.
* `A2` Effect allele (encoded as 1 in PLINK bed file).
* `MAF` Minor allele frequency
* `NCHROBS` Number of samples times 2.

Minor allele frequency files for the East Asian and European populations can
be downloaded [here](https://data.broadinstitute.org/alkesgroup/S-LDXR/maf.tar.gz).
