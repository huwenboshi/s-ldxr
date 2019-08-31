# Estimating stratified squared trans-ethnic genetic correlation and its
enrichment

This page describes how to estimate stratified squared trans-ethnic genetic
correlation, \\(r^2_{g}(C)\\), and its enrichment, \\(\lambda^2(C)\\).

### The command

X-LDSC estimates \\(r^2_{g}(C)\\) and \\(\lambda^2(C)\\) with the
following command.

```
python <software directory>/x-ldsc.py \
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

### The output

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
