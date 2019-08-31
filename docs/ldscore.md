# Estimating LD scores

This page describes how to estimate LD scores

### Running the tool

X-LDSC estimates LD scores with the following command.

```
for chrom in $(seq 22)
do
    python <software directory>/x-ldsc.py \
        --score allelic \
        --ld-wind-cm 1.0 \
        --print-snps <a list of SNPs to print> \
        --bfile <EAS reference panel directory>/1000G.EAS.${chrom} \
                <EUR reference panel directory>/1000G.EUR.${chrom} \
        --annot <annotation directory>/baseline.${chrom}.annot.gz \
        --out <output directory>/EAS_EUR_baseline_chr${chrom}
done
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: The above for loop can be parallelized by chromosome, if a
computer cluster is available.
</div>

Here are the meanings of the flags:

* `--score` specifies the type of LD scores to be estimated. Here, one should
almost always put "allelic" as the score type for estimating per-allele effect
correlation. Although "standardized" score type is also supported, we do not
recommend to only use it in computing regression weights, as trans-ethnic
genetic correlation of standardized causal effect sizes is not very
interpretable.

* `--ld-wind-cm` specifies the maximum window size in centimorgan for
estimating LD scores. The default and recommended value is 1.0.

* `--print-snps` specifies the file that contains a list of SNPs in plain text
without header, for which the LD scores are to be printed. We recommend to
print all SNPs with minor allele frequency greater than 1% in both populations.

* `--bfile` takes two argument -- reference panel for population 1, and
reference panel for population 2. All reference panels should be in PLINK
format, and have the same set of SNPs.

* `--annot` specifies the annotation file.

* `--out` specifies prefix of the output files.

### Interpreting the output

After executing the command above, 4 files will be created for each
chromosome (i.e. 88 files for all 22 chromosomes):

* **step1_chr22.info.gz** - contains the locus information, including start and end
positions, number of SNPs, rank of LD matrices, and sample size

* **step1_chr22.log** - contains helpful information for debugging, including
number of SNPs, number of SNPs filtered, etc.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: Log files are very useful in pinpointing bugs of the
software. Please include the log file in the email in any bug report.
</div>
