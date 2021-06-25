# Estimating enrichment of stratified squared trans-ethnic genetic correlation

This page describes the steps to estimate enrichment of stratified squared
trans-ethnic genetic correlation,
\\(\lambda^2(C) ={ {r^2_g(C)} \over {r^2_g} }\\), the ratio between squared
trans-ethnic genetic correlation of annotation \\( C \\) and genome-wide
squared trans-ethnic genetic correlation.

## Typical command

S-LDXR estimates \\(\lambda^2(C)\\) with the following command.

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

This command typically takes 10 to 15 minutes to run on a stand alone computer.

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

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: By default, S-LDXR corrects for bias in ratio estimation
analytically. We also provide the option to correct for bias using jackknife.
This can be achieved by adding the --use-jackknife-bias-adj flag.
</div>

## Output

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

## Estimating \\( \lambda^2(C) \\) for continuous-valued annotations

The following command estimates enrichment of stratified squared trans-ethnic
genetic correlation for quintiles of continuous-valued annotations.

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

This step typically takes 2 to 5 minutes to run on a stand alone computer.

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

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: By default, S-LDXR corrects for bias in ratio estimation
analytically. We also provide the option to correct for bias using jackknife.
This can be achieved by adding the --use-jackknife-bias-adj flag.
</div>

## Expected \\(\lambda^2(C)\\) from continuous-valued annotations

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

This command typically takes 2 to 5 minutes to run on a stand alone computer.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: It is important to include the base (not baseline) annotation.
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

This command typically takes 2 to 5 minutes to run on a stand alone computer.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: It is important to include the base (not baseline) annotation.
By default, S-LDXR corrects for bias in ratio estimation
analytically. We also provide the option to correct for bias using jackknife.
This can be achieved by adding the --use-jackknife-bias-adj flag.
</div>

The output is the same as that of the command for continuous-valued
annotations.

## Interpreting the output

The output files of S-LDXR contain the following columns.

1. `ANNOT` name of the annotation

2. `NSNP` number of SNPs for binary annotations (sum of annotation values for
continuous-valued annotations)
3. `STD` standard deviation of the annotation across SNPs

4. `TAU1` heritability annotation coefficient of population 1

5. `TAU1_SE` standard error heritability annotation coefficient of population 1

6. `TAU2` heritability annotation coefficient of population 2

7. `TAU2_SE` standard error heritability annotation coefficient of population 2

8. `THETA` trans-ethnic genetic covariance annotation coefficient

9. `THETA_SE` standard error of trans-ethnic genetic covariance annotation coefficient

10. `HSQ1` stratified heritability in population 1

11. `HSQ1_SE` standard error of stratified heritability in population 1

12. `HSQ2` stratified heritability in population 2

13. `HSQ2_SE` standard error of stratified heritability in population 2

14. `GCOV` stratified trans-ethnic genetic covariance

15. `GCOV_SE` standard error of stratified trans-ethnic genetic covariance

16. `GCOR` stratified trans-ethnic genetic correlation

17. `GCOR_SE` standard error for the estimated stratified trans-ethnic
genetic correlation

18. `GCORSQ` stratified squared trans-ethnic genetic correlation

19. `GCORSQ_SE` standard error of stratified squared trans-ethnic genetic
correlation 

20. `HSQ1_ENRICHMENT` heritability enrichment in population 1

21. `HSQ1_ENRICHMENT_SE` standard error of heritability enrichment in
population 1

22. `HSQ2_ENRICHMENT` heritability enrichment in population 2

23. `HSQ2_ENRICHMENT_SE` standard error of heritability enrichment in
population 2

24. `GCOV_ENRICHMENT` genetic covariance enrichment

25. `GCOV_ENRICHMENT_SE` standard error of genetic covariance enrichment

26. `GCORSQ_ENRICHMENT` estimated enrichment of stratified squared trans-ethnic
genetic correlation enrichment

27. `GCORSQ_ENRICHMENT_SE` standard error of estimated enrichment of stratified
squared trans-ethnic genetic correlation

28. `GCORSQ_ENRICHMENT_P` p-value for testing whether enrichment of stratified
trans-ethnic genetic correlation is different from 1. Here the p-value
is obtained from a t distribution with degree of freedom equal to the number
of jackknife blocks minus one, where the test statistic is
\\( { {\hat{\lambda}^2(C)} \over {s.e.(\hat{\lambda}^2(C)) } }\\).

29. `GCOVSQ_DIFF` estimated
\\( \hat{D}^2(C) = \hat{\rho}^2_g(C) - \hat{r}^2_g \hat{h}^2_g(C) \hat{h}^2_g(C) \\),
the difference between stratified squared trans-ethnic genetic covariance
of annotation \\( C \\), and \\( \hat{r}^2_g \hat{h}^2_g(C) \hat{h}^2_g(C) \\),
the expected squared trans-ethnic genetic covariance based on genome-wide
squared trans-ethnic genetic correlation and heritabilities.

30. `GCOVSQ_DIFF_SE` standard error for the estimated \\( \hat{D}^2(C) \\)

31. `GCOVSQ_DIFF_P` p-value for testing whether \\( \hat{D}^2(C) \\) is different
from 0, obtained from a t distribution with degree of freedom equal to the number
of jackknife blocks minus one, where the test statistic is
\\( { {\hat{D}^2(C)} \over {s.e.(\hat{D}^2(C)) } }\\). This test is equivalent
to testing whether \\( \lambda^2(C) \\) is different from 1. But the p-value
is better calibrated.

<div style="background-color:rgba(230, 230, 250, 1.0);">
<b>Note</b>: We recommend to use the column GCORSQ_ENRICHMENT and
GCORSQ_ENRICHMENT_SE for meta-analysis of results across traits, and use
GCOVSQ_DIFF_P to test for enrichment/depletion.
</div>

