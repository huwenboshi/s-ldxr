# Estimating stratified squared trans-ethnic genetic correlation and its enrichment

This page describes how to estimate stratified squared trans-ethnic genetic
correlation, \\(r^2_{g}(C)\\), and its enrichment, \\(\lambda^2(C)\\).

### Typical command

S-LDXR estimates \\(r^2_{g}(C)\\) and \\(\lambda^2(C)\\) with the
following command.

```
python <software directory>/s-ldxr.py \
    --gcor <summary stats directory for EAS>/EAS_sumstats.gz \
           <summary stats directory for EUR>/EUR_sumstats.gz \
    --ref-ld-chr <baseline LD score directory>/EAS_EUR_baseline_chr \
                 <AVGLLD LD score directory>/EAS_EUR_avglld_chr \
                 <BSTAT LD score directory>/EAS_EUR_bstat_chr \
                 <ALLELEAGE LD score directory>/EAS_EUR_alleleage_chr \
    --w-ld-chr <regression weight directory>/EAS_EUR_weight_chr \
    --frqfile <EAS MAF directory>/1000G.EAS. \
              <EUR MAF directory>/1000G.EUR. \
    --annot <baseline annotation directory>/baseline. \
            <AVGLLD annotation directory>/avglld. \
            <BSTAT annotation directory>/bstat. \
            <ALLELEAGE annotation directory>/alleleage. \
    --apply-shrinkage 0.5 \
    --save-pseudo-coef \
    --out TRAIT_EAS_EUR.txt
```

Here are the meanings of the flags:

* `--gcor` specifies the summary stats files. This flag takes 2 arguments -
summary stats for population 1 and summary stats for population 2.

* `--ref-ld-chr` specifies prefix of the LD score files. This flag takes one
or more arguments – one may put as many LD score files as one wishes.

* `--w-ld-chr` specifies prefix of the regression weights. These are
standardized LD scores calculated from regression SNPs.

* `--frqfile` specifies prefix of minor allele frequency files.

* `--annot` specifies prefix of the annotation files. This flags also takes
one or more arguments.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: The order one specifies the annotation files must
be the same as the order one specifies the LD score files. The annotation
files must also be the same files that one uses to obtain the LD scores.
</div>

* `--apply-shrinkage` adjusts the level of shrinkage (the \\(\alpha\\) tuning
parameter in the paper). This should be a number between 0 and 1.

* `--save-pseudo-coef` If this flag is specified, jackknife pseudo values of
the coefficients will be saved. This flag is optional.

* `--out` specifies the output file name.

### Output

After executing the above command, 5 files will be generated.

* `TRAIT_EAS_EUR.txt` output file containing the estimates.

* `TRAIT_EAS_EUR.txt.log` log file containing information for debugging.

* `TRAIT_EAS_EUR.txt.pseudo_tau1.gz` jackknife pseudo values for \\(\tau_C\\)
coefficients for population 1.

* `TRAIT_EAS_EUR.txt.pseudo_tau2.gz` jackknife pseudo values for \\(\tau_C\\)
coefficients for population 2.

* `TRAIT_EAS_EUR.txt.pseudo_theta.gz` jackknife pseudo values for
\\(\theta_C\\) genetic covariance coefficients.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: The pseudo coefficients will only be saved if the
--save-pseudo-coef flag is specified.
</div>

### Continuous-valued annotations

The following command estimates enrichment of stratified squared trans-ethnic
genetic correlation of quintiles of continuous-valued annotations.

```
python <software directory>/cont_annot_gcor.py \
    --coef TRAIT_EAS_EUR.txt \
    --frqfile <EAS MAF directory>/1000G.EAS. \
              <EUR MAF directory>/1000G.EUR. \
    --annot <baseline annotation directory>/baseline. \
            <AVGLLD annotation directory>/avglld. \
            <BSTAT annotation directory>/bstat. \
            <ALLELEAGE annotation directory>/alleleage. \
    --names AVGLLD BSTAT ALLELEAGE \
    --nbins 5 \
    --out TRAIT_EAS_EUR_contannot.txt
```

Here are the meanings of the flags.

* `--coef` specifies the output from the previous step. The jackknife pseudo
coefficients will be loaded automatically.

* `--frqfile` specifies prefix of minor allele frequency files.

* `--annot` specifies prefix of the annotation files. This flags also takes
one or more arguments.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: The order one specifies the annotation files must be the same as the
order of annotations in TRAIT_EAS_EUR.txt.
</div>

* `--names` specifies the names of the continuous annotations for which one
wishes to compute enrichment at quintiles.

* `--nbins` specifies the number of bins to bin the SNPs based on the values
of their continuous annotation. The default is 5 (i.e. quintiles).

* `--out` specifies the output file name.

Additionally, users may use the `--apply-shrinkage` flag to adjust the level
of shrinkage.

After executing the above command, 2 files will be created.

* `TRAIT_EAS_EUR_contannot.txt` contains the estimates.

* `TRAIT_EAS_EUR_contannot.txt.log` is the log file for debugging purpose.

### Expected \\(r^2_g(C)\\) and \\(\lambda^2(C)\\) from continuous-valued annotations

Estimating expected \\(r^2_g(C)\\) and \\(\lambda^2(C)\\) from
continuous-valued annotations requires two steps.

The first step gets the coefficients (\\(\tau_{1C}\\), \\(\tau_{2C}\\),
and \\(\theta_{C}\\)) of each continuous-valued annotations

```
python <software directory>/s-ldxr.py \
    --gcor <summary stats directory for EAS>/EAS_sumstats.gz \
           <summary stats directory for EUR>/EUR_sumstats.gz \
    --ref-ld-chr <base LD score directory>/EAS_EUR_allelic_chr \
                 <AVGLLD LD score directory>/EAS_EUR_allelic_chr \
                 <BSTAT LD score directory>/EAS_EUR_allelic_chr \
                 <ALLELEAGE LD score directory>/EAS_EUR_allelic_chr \
    --w-ld-chr <regression weight directory>/EAS_EUR_weight_chr \
    --frqfile <EAS MAF directory>/1000G.EAS. \
              <EUR MAF directory>/1000G.EUR. \
    --annot <base annotation directory>/base. \
            <AVGLLD annotation directory>/avglld. \
            <BSTAT annotation directory>/bstat. \
            <ALLELEAGE annotation directory>/alleleage. \
    --save-pseudo-coef \
    --out ./TRAIT_EAS_EUR_step1.txt
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: It is important to always include the base (not
baseline) annotation.
</div>

The output is the same as that of a typical command.

The second step uses the coefficients from the first step to obtain expected
\\(r^2_g(C)\\) and \\(\lambda^2(C)\\) from continuous-valued annotations.

```
python <software directory>/pred_binannot_from_contannot.py \
    --coef ./TRAIT_EAS_EUR_step1.txt \
    --frqfile <EAS MAF directory>/1000G.EAS. \
              <EUR MAF directory>/1000G.EUR. \
    --cont-annot <base annotation directory>/base. \
                 <AVGLLD annotation directory>/avglld. \
                 <BSTAT annotation directory>/bstat. \
                 <ALLELEAGE annotation directory>/alleleage. \
    --bin-annot <base annotation directory>/base. \
                <binary annotation directory>/annot_name. \
    --apply-shrinkage 0.5 \
    --out ./TRAIT_EAS_EUR_step2.txt
```

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: It is important to always include the base (not
baseline) annotation.
</div>

The output is the same as that of the command for continuous-valued
annotations.

### Interpreting the output

The output files of S-LDXR contain the following columns.

1. `ANNOT` name of the annotation

2. `NSNP` number of SNPs for binary annotations (sum of annotation values for
continuous-valued annotations)
3. `STD` standard deviation of the annotation across SNPs

4. `TAU1` heritability coefficient of population 1

5. `TAU1_SE` standard error heritability coefficient of population 1

6. `TAU2` heritability coefficient of population 2

7. `TAU2_SE` standard error heritability coefficient of population 2

8. `THETA` trans-ethnic genetic covariance coefficient

9. `THETA_SE` standard error of trans-ethnic genetic covariance coefficient

10. `HSQ1` heritability in population 1

11. `HSQ1_SE` standard error of heritability in population 1

12. `HSQ2` heritability in population 2

13. `HSQ2_SE` standard error of heritability in population 2

14. `GCOV` trans-ethnic genetic covariance

15. `GCOV_SE` standard error of trans-ethnic genetic covariance

16. `GCORSQ` stratified squared trans-ethnic genetic correlation

17. `GCORSQ_SE` standard error of stratified squared trans-ethnic genetic
correlation 

18. `HSQ1_ENRICHMENT` heritability enrichment in population 1

19. `HSQ1_ENRICHMENT_SE` standard error of heritability enrichment in
population 1

20. `HSQ2_ENRICHMENT` heritability enrichment in population 2

21. `HSQ2_ENRICHMENT_SE` standard error of heritability enrichment in
population 2

22. `GCOV_ENRICHMENT` genetic covariance enrichment

23. `GCOV_ENRICHMENT_SE` standard error of genetic covariance enrichment

24. `GCORSQ_ENRICHMENT` estimated enrichment of stratified squared trans-ethnic
genetic correlation enrichment

25. `GCORSQ_ENRICHMENT_SE` standard error of estimated enrichment of stratified
squared trans-ethnic genetic correlation
