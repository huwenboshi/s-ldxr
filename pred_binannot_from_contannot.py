#!/usr/bin/python
# (c) 2018-2023 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from src.estimation import *

# main function
def main():

    # get command line argument and initialize log
    args = get_command_line()

    init_log(args)

    # parse command line 
    argmap = check_command_line(args)

    # load coefficients
    coef = pd.read_table(argmap['coef'])
    tau1 = coef['TAU1'].values.astype(np.float32)
    tau2 = coef['TAU2'].values.astype(np.float32)
    theta = coef['THETA'].values.astype(np.float32)

    # load pseudo coefficients
    ps_tau1 = load_pseudo_coef(argmap['coef']+'.pseudo_tau1.gz')
    ps_tau2 = load_pseudo_coef(argmap['coef']+'.pseudo_tau2.gz')
    ps_theta = load_pseudo_coef(argmap['coef']+'.pseudo_theta.gz')
    nblock = ps_tau1.shape[0]

    # load ld-related annotations
    start_chrom, stop_chrom = argmap['use-chrom']
    annot_snps, all_ld_annot, ld_annot_mat = load_annot(argmap['ld-annot'],
        start_chrom, stop_chrom)

    # load baseline annotations
    _, all_bl_annot, bl_annot_mat = load_annot(argmap['baseline-annot'],
        start_chrom, stop_chrom)

    # load frqfile
    freq1 = load_frqfile(argmap['frqfile'][0], start_chrom, stop_chrom)
    freq2 = load_frqfile(argmap['frqfile'][1], start_chrom, stop_chrom)

    # use frq to filter out snps for estimation
    min_maf = argmap['min-maf']
    annot_idx = intersect_snp_estimation(annot_snps['SNP'], freq1,
        freq2, min_maf)
    ld_annot_mat = ld_annot_mat[annot_idx,:] 
    bl_annot_mat = bl_annot_mat[annot_idx,:]

    # get summary information for annotation matrix
    annot_nsnp, bl_annot_std = get_annot_sumstat(bl_annot_mat)

    # get hsq and gcov for each annotation
    hsq1,hsq1_se,ps_hsq1 = get_sum_bin(ld_annot_mat,tau1,ps_tau1,bl_annot_mat)
    hsq2,hsq2_se,ps_hsq2 = get_sum_bin(ld_annot_mat,tau2,ps_tau2,bl_annot_mat)
    gcov,gcov_se,ps_gcov=get_sum_bin(ld_annot_mat,theta,ps_theta,bl_annot_mat)
    hsq_jkcov = get_jkcov(hsq1, ps_hsq1, hsq2, ps_hsq2)

    # get enrichment estimate
    hsq1_en,hsq1_en_se,ps_hsq1_en = get_enrichment(hsq1, ps_hsq1, annot_nsnp)
    hsq2_en,hsq2_en_se,ps_hsq2_en = get_enrichment(hsq2, ps_hsq2, annot_nsnp)
    gcov_en,gcov_en_se,ps_gcov_en = get_enrichment(gcov, ps_gcov, annot_nsnp)

    # get sqaured genetic correlation estimates
    gcorsq, gcorsq_se, ps_gcorsq = get_gcorsq(hsq1, ps_hsq1, hsq2, ps_hsq2,
        gcov, ps_gcov, annot_nsnp, False, argmap['apply-shrinkage'])

    # estimate enrichment of squared genetic correlation
    gcorsq_en, gcorsq_en_se = get_gcorsq_enrichment(gcorsq, ps_gcorsq)
    tstat = np.fabs((gcorsq_en-1.0) / (gcorsq_en_se+1e-16))
    gcorsq_en_pval = (1.0-scipy.stats.t.cdf(tstat, nblock-1))*2.0

    # difference in gcov squared
    gcovsq_diff, gcovsq_diff_se, ps_gcovsq_diff = get_gcovsq_diff(hsq1,
        ps_hsq1, hsq2, ps_hsq2, gcov, gcov_se, ps_gcov, gcorsq, ps_gcorsq)
    tstat = np.fabs(gcovsq_diff / (gcovsq_diff_se+1e-16))
    gcovsq_diff_pval = (1.0-scipy.stats.t.cdf(tstat, nblock-1))*2.0

    # collect results in a data frame
    out_df = pd.DataFrame()
    out_df['ANNOT'] = all_bl_annot['ANNOT']
    out_df['NSNP'] = annot_nsnp.astype(int)
    out_df['STD'] = bl_annot_std
    out_df['HSQ1'] = hsq1; out_df['HSQ1_SE'] = hsq1_se
    out_df['HSQ2'] = hsq2; out_df['HSQ2_SE'] = hsq2_se
    out_df['GCOV'] = gcov; out_df['GCOV_SE'] = gcov_se
    out_df['GCORSQ'] = gcorsq; out_df['GCORSQ_SE'] = gcorsq_se
    out_df['HSQ1_ENRICHMENT'] = hsq1_en
    out_df['HSQ1_ENRICHMENT_SE'] = hsq1_en_se
    out_df['HSQ2_ENRICHMENT'] = hsq2_en
    out_df['HSQ2_ENRICHMENT_SE'] = hsq2_en_se
    out_df['GCOV_ENRICHMENT'] = gcov_en
    out_df['GCOV_ENRICHMENT_SE'] = gcov_en_se
    out_df['GCORSQ_ENRICHMENT'] = gcorsq_en
    out_df['GCORSQ_ENRICHMENT_SE'] = gcorsq_en_se
    out_df['GCORSQ_ENRICHMENT_P'] = gcorsq_en_pval
    out_df['GCOVSQ_DIFF'] = gcovsq_diff
    out_df['GCOVSQ_DIFF_SE'] = gcovsq_diff_se
    out_df['GCOVSQ_DIFF_P'] = gcovsq_diff_pval
    
    # write results to file
    out_df.to_csv(argmap['out'], sep='\t', index=False)

    end_log()

def get_sum_bin(annot_mat, coef, ps_coef, bin_annot_mat):
    """
    Estimate the heritability / genetic covariance in each functional
    annotation
    """

    # get per snp variance / covariance estimate
    nannot = bin_annot_mat.shape[1]
   
    # exclude coef for the intercept 
    annot_cov = np.dot(bin_annot_mat.T, annot_mat)    
    annot_est = np.dot(annot_cov, coef)

    # get standard error
    nblock = ps_coef.shape[0]
    all_ps_annot_est = np.zeros((nblock, nannot), dtype=np.float32)
    for j in range(nblock):
        all_ps_annot_est[j,:] = np.dot(annot_cov, ps_coef[j][:-1])

    # get se
    nblock = np.float32(nblock)
    diffsq = np.square(all_ps_annot_est-annot_est)
    annot_est_se = np.sqrt((nblock-1)*np.mean(diffsq,axis=0))

    return annot_est, annot_est_se, all_ps_annot_est

# load jack knife coefficient
def load_pseudo_coef(filename):
    ps_coef = np.loadtxt(filename, dtype=np.float32)
    return ps_coef

# initialize log
def init_log(args):

    # get log file name
    log_file_name = args.out
    log_file_name += '.log'

    # create the log file
    log_format = '[%(levelname)s] %(message)s'
    logging.basicConfig(filename=log_file_name, filemode="w",
        level=logging.DEBUG, format=log_format)

    # add stderr as a stream handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(log_format)
    stderr_handler.setFormatter(formatter)
    logging.getLogger().addHandler(stderr_handler)

    # log time and command issued
    #logging.info(title_str)
    specified = set([val for val in sys.argv if val[0] == '-'])
    cmd_str = sys.argv[0] + ' \\\n'
    for arg in vars(args):
        if ('--' + arg) in specified:
            param = getattr(args, arg)
            if type(param) == list:
                param = ' '.join([str(p) for p in param])
            elif type(param) == bool:
                param = ''
            cmd_str += '        --{} {} \\\n'.format(arg, param)
    cmd_str = cmd_str.strip()[0:len(cmd_str)-3]
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command started at: %s' % cur_time)
    logging.info('Command issued:\n    {}'.format(cmd_str))

# end the log
def end_log():
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command finished at: %s' % cur_time)

# check command line
def check_command_line(args):

    # parse command line
    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    return argmap

# get command line input
def get_command_line():
    
    # create the help document
    parser = argparse.ArgumentParser(description='Partition trans-ethnic '
        'genetic correlation by functional categories')

    parser.add_argument('--coef', dest='coef', type=str, required=False,
        help='Coefficient files')

    parser.add_argument('--ld-annot', dest='ld-annot', type=str, required=False,
        nargs='+', help='LD-related annotation file')

    parser.add_argument('--baseline-annot', dest='baseline-annot', type=str,
        required=False, nargs='+', help='Baseline annotation file')

    parser.add_argument('--use-chrom', dest='use-chrom', type=int,
        required=False, default=[1, 22], nargs=2,
        help='Specific which chromosomes to use')
   
    parser.add_argument('--frqfile', dest='frqfile', type=str, nargs=2,
        required=False, help='Prefix of the frequency files')

    parser.add_argument('--min-maf', dest='min-maf', type=float, default=0.05,
        required=False, help='Minimum MAF for performing the estimation')

    parser.add_argument('--apply-shrinkage', dest='apply-shrinkage',
        type=float, default=0.5, required=False,
        help='Shrinkage applied when estimating genetic correlation')

    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
