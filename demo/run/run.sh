frqdir=/n/groups/price/huwenbo/PART_TE_GCOV/simulation/data/haplotype/plink/relcutoff_35k
annot_dir=/n/groups/price/huwenbo/PART_TE_GCOV/simulation/data/annotation_35k

python ../../s-ldxr.py \
    --gcor ../data/sim1_pop1.sumstats.gz \
           ../data/sim1_pop2.sumstats.gz \
    --ref-ld-chr $annot_dir/annot_baseline/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_gerp/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_fst/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_recomb/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_diversity/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_avglld/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_bstat/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_alleleage/ldscore_allelic/EAS_EUR_allelic_chr \
                 $annot_dir/annot_cpg/ldscore_allelic/EAS_EUR_allelic_chr \
    --w-ld-chr $annot_dir/annot_weight_5pct/ldscore_standardized/EAS_EUR_standardized_chr \
    --frqfile $frqdir/EAS_ref/EAS_ref. $frqdir/EUR_ref/EUR_ref. \
    --annot $annot_dir/annot_baseline/annot_files/baseline. \
            $annot_dir/annot_gerp/annot_files/gerp. \
            $annot_dir/annot_fst/annot_files/fst. \
            $annot_dir/annot_recomb/annot_files/recomb. \
            $annot_dir/annot_diversity/annot_files/diversity. \
            $annot_dir/annot_avglld/annot_files/avglld. \
            $annot_dir/annot_bstat/annot_files/bstat. \
            $annot_dir/annot_alleleage/annot_files/alleleage. \
            $annot_dir/annot_cpg/annot_files/cpg. \
    --use-chrom 1 3 \
    --n-blocks 50 \
    --min-maf 0.05 \
    --apply-shrinkage 0.5 \
    --save-pseudo-coef \
    --out ./sim1.txt

python ../../cont_annot_gcor.py \
    --coef ./sim1.txt \
    --frqfile ${frqdir}/EAS_ref/EAS_ref. ${frqdir}/EUR_ref/EUR_ref. \
    --annot $annot_dir/annot_baseline/annot_files/baseline. \
            $annot_dir/annot_gerp/annot_files/gerp. \
            $annot_dir/annot_fst/annot_files/fst. \
            $annot_dir/annot_recomb/annot_files/recomb. \
            $annot_dir/annot_diversity/annot_files/diversity. \
            $annot_dir/annot_avglld/annot_files/avglld. \
            $annot_dir/annot_bstat/annot_files/bstat. \
            $annot_dir/annot_alleleage/annot_files/alleleage. \
            $annot_dir/annot_cpg/annot_files/cpg. \
    --names GERP_NS FST Recomb_rate_10kb Nucleotide_Diversity_10kb AVGLLD BSTAT ALLELEAGE CpG \
    --use-chrom 1 3 \
    --nbins 5 \
    --min-maf 0.05 \
    --apply-shrinkage 0.5 \
    --out ./sim1_cont5.txt
