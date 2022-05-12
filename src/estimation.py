import sys, logging
import numpy as np
import pandas as pd
import scipy, scipy.linalg
import statsmodels.api as sm
from .sumstats import *
from .score import *
from .annot import *

def intersect_snp_regression(snp1, snp2, score_snps, regsnps, frqsnps):
    """
    Intersect SNPs for gression
    """

    # get shared snps
    snp12 = snp1[snp1.isin(snp2)]
    snp12_score = snp12[snp12.isin(score_snps)]
    shared = snp12_score[snp12_score.isin(regsnps)]

    # get index
    snp1_idx = np.where(snp1.isin(shared))[0]
    snp2_idx = np.where(snp2.isin(shared))[0]
    score_idx = np.where(score_snps.isin(shared))[0]
    reg_idx = np.where(regsnps.isin(shared))[0]
    freq_idx = np.where(frqsnps.isin(shared))[0]

    return snp1_idx, snp2_idx, score_idx, reg_idx, freq_idx

def intersect_snp_estimation(annot_snps, freq1, freq2, min_maf):
    """
    Intersect SNPS for estimation
    """
    annot_idx = np.where((freq1['MAF']>min_maf) & 
                         (freq2['MAF']>min_maf))[0]
    return annot_idx

def create_block(start_idx, stop_idx, nblock):
    """
    Create blocks for block jackknife
    """

    block_size = int(np.ceil(float(stop_idx-start_idx+1)/float(nblock)))
    cuts = range(start_idx, stop_idx, block_size)
    blocks = []
    
    for i in range(len(cuts)-1):
        start = cuts[i]
        stop = cuts[i+1]
        blocks.append(np.arange(start, stop))
        
    blocks.append(np.arange(cuts[len(cuts)-1], stop_idx+1))

    return blocks

def regression(x, y, nblock, fit_intercept):
    """
    Perform least square regression with block jack knife
    """

    # perform regression
    xtx = np.dot(x.T, x)
    xty = np.dot(x.T, y)
    ncoef = xtx.shape[0]

    # obtain coefficient
    coef = np.zeros(ncoef, dtype=np.float32)
    if fit_intercept == 'no':
        coef[:-1] = np.linalg.solve(xtx[:-1,:-1], xty[:-1])
    else:
        coef = np.linalg.solve(xtx, xty)

    if nblock == 0:
        return coef

    # get jacknife estimate
    all_blocks = create_block(0, x.shape[0]-1, nblock)
    ps_coef = np.zeros((len(all_blocks), x.shape[1]), dtype=np.float32)
    for i in range(len(all_blocks)):
        block = all_blocks[i]
        xtx_block = xtx - np.dot(x[block,:].T, x[block,:])
        xty_block = xty - np.dot(x[block,:].T, y[block])
        coef_block = np.zeros(ncoef, dtype=np.float32)
        if fit_intercept == 'no':
            coef_block[:-1]=np.linalg.solve(xtx_block[:-1,:-1],xty_block[:-1])
        else:
            coef_block = np.linalg.solve(xtx_block, xty_block)
        ps_coef[i,:] = coef_block

    # get se
    nblock = np.float32(nblock)
    mean_ps_coef = np.mean(ps_coef, axis=0)
    se = np.sqrt((nblock-1)*np.mean(np.square(ps_coef-mean_ps_coef), axis=0))

    return coef, se, ps_coef

def robust_regression(x, y, nblock, fit_intercept):
    """
    Perform least square regression with block jack knife
    """

    # tuning constant
    cval = 10.0

    # perform regression
    ncoef = x.shape[1]

    # obtain coefficient
    coef = np.zeros(ncoef, dtype=np.float32)
    if fit_intercept == 'no':
        model = sm.RLM(y, x[:,:-1], M=sm.robust.norms.TukeyBiweight(c=cval))
        results = model.fit()
        coef[:-1] = results.params
    else:
        model = sm.RLM(y, x, M=sm.robust.norms.TukeyBiweight(c=cval))
        results = model.fit()
        coef = results.params

    if nblock == 0:
        return coef

    # get jacknife estimate
    all_blocks = create_block(0, x.shape[0]-1, nblock)
    all_idx = np.arange(0, x.shape[0])
    ps_coef = np.zeros((len(all_blocks), x.shape[1]), dtype=np.float32)
    for i in range(len(all_blocks)):
        block = all_blocks[i]
        use_idx = np.delete(all_idx, block)
        x_block = x[use_idx,:]
        y_block = y[use_idx]
        coef_block = np.zeros(ncoef, dtype=np.float32)
        
        if fit_intercept == 'no':
            model = sm.RLM(y_block, x_block[:,:-1],
                M=sm.robust.norms.TukeyBiweight(c=cval))
            results = model.fit()
            coef_block[:-1] = results.params
        else:
            model = sm.RLM(y_block,x_block,
                M=sm.robust.norms.TukeyBiweight(c=cval))
            results = model.fit()
            coef_block = results.params
        
        ps_coef[i,:] = coef_block

    # compute standard error
    mean_ps_coef = np.mean(ps_coef, axis=0)
    se = np.sum(np.square(ps_coef-mean_ps_coef), axis=0)
    se = np.sqrt(se*np.float32(nblock-1)/np.float32(nblock))

    return coef, se, ps_coef

def get_coef(score, zsc1, zsc2, w, nblock,
    fit_intercept, use_robust_regression, subtract):
    """
    Obtain coefficient for heritability / genetic covariance
    """
    
    # set up the regression
    nsnp = score.shape[0]
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    if fit_intercept == 'no':
        zprod -= subtract
    score[:,:-1] = score[:,:-1] * nprod[:,np.newaxis]
    
    # scale the matrix to improve matrix condition
    nbar = np.mean(nprod)
    score[:,:-1] /= nbar
    
    # apply the weight
    score_w = score * np.sqrt(w[:,np.newaxis])
    zprod_w = zprod * np.sqrt(w)

    # run regression with jackknife
    if nblock == 0:
        if use_robust_regression == False:
            coef = regression(score_w, zprod_w, 0, fit_intercept)
        else:
            coef = robust_regression(score, zprod, 0, fit_intercept)
        coef[:-1] /= nbar
        return coef
    
    if use_robust_regression == False:
        coef, coef_se, ps_coef = regression(score_w, zprod_w,
            nblock, fit_intercept)
    else:
        coef, coef_se, ps_coef = robust_regression(score, zprod,
            nblock, fit_intercept)
        
    # rescale the coefficient
    coef[:-1] /= nbar
    coef_se[:-1] /= nbar
    ps_coef[:,:-1] /= nbar

    return np.float64(coef), np.float64(coef_se), np.float64(ps_coef)

def get_coef_raw(ldscore_mat, zsc1, zsc2, max_int):
    """
    Obtain coefficient for heritability / genetic covariance
    """
   
    # get dimension
    nsnp,ncoef = ldscore_mat.shape[0],ldscore_mat.shape[1]

    # set up variables
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    score_sum = np.sum(ldscore_mat[:,:-1], axis=1)
    
    # get averages
    mean_ldscore = np.mean(score_sum)
    mean_zprod = np.mean(zprod)
    mean_nprod = np.mean(nprod)

    # estimate intercept
    zprod_sorted = np.sort(zprod)
    idx = int(float(nsnp)*0.95)
    intercept = np.mean(zprod_sorted[0:idx])
    if intercept > max_int:
        intercept = max_int

    # get raw coef
    coef = (mean_zprod-intercept) / (mean_nprod*mean_ldscore)

    return coef, intercept

# get the sum of variance / covariance for each annotation
def get_sum(annot_mat, coef, ps_coef):
    """
    Estimate the heritability / genetic covariance in each functional
    annotation
    """

    # get per snp variance / covariance estimate
    nannot = annot_mat.shape[1]
   
    # exclude coef for the intercept 
    annot_cov = np.dot(annot_mat.T, annot_mat)    
    annot_est = np.dot(annot_cov, coef[:-1])

    annot_cov = annot_cov.astype(np.float64)
    annot_est = annot_est.astype(np.float64)
    
    # get standard error
    nblock = ps_coef.shape[0]
    all_ps_annot_est = np.zeros((nblock, nannot), dtype=np.float64)
    for j in range(nblock):
        all_ps_annot_est[j,:] = np.dot(annot_cov, ps_coef[j][:-1])
   
    # get se
    nblock = np.float64(nblock)
    mean_all_ps_annot_est = np.mean(all_ps_annot_est, axis=0)
    diffsq = np.square(all_ps_annot_est-mean_all_ps_annot_est)
    annot_est_se = np.sqrt((nblock-1)*np.mean(diffsq,axis=0))

    return annot_est.astype(np.float64), annot_est_se.astype(np.float64), \
        all_ps_annot_est.astype(np.float64)

# get jackknife covariance of heritability estimates
def get_jkcov(est1, ps_est1, est2, ps_est2):
    """
    Use jack knife pseudo values of heritability estimation to obtain
    jack knife covariance
    """

    nblock = ps_est1.shape[0]
    mean_ps_est1 = np.mean(ps_est1, axis=0)
    mean_ps_est2 = np.mean(ps_est2, axis=0)
    jkcov = np.mean((ps_est1-mean_ps_est1)*(ps_est2-mean_ps_est2), axis=0)
    jkcov *= (np.float64(nblock)-1.0)
    jkcov = jkcov.astype(np.float64)

    return jkcov

def get_enrichment(est, ps_est, annot_nsnp):
    
    """
    Estimate enrichment of functional annotation
    """

    # get dimension
    nannot = annot_nsnp.shape[0]
    tot_nsnp = np.float64(annot_nsnp[0])
    annot_en = est*tot_nsnp / (est[0] * annot_nsnp)

    # get standard error
    nblock = ps_est.shape[0]
    all_ps_annot_en = np.zeros((nblock, nannot), dtype=np.float64)
    for j in range(nblock):
        all_ps_annot_en[j,:] = ps_est[j,:]*tot_nsnp/(ps_est[j][0]*annot_nsnp)

    # adjust for bias using jackknife
    jknife_mean = np.mean(all_ps_annot_en, axis=0)
    bias = (nblock-1) * (annot_en - jknife_mean)
    annot_en_adj = annot_en - bias
    factor = np.float64(nblock)/(np.float64(nblock)-1.0)
    all_ps_annot_en_adj = all_ps_annot_en - factor*bias

    # get se
    nblock = np.float64(nblock)
    mean_all_ps_annot_en_adj = np.mean(all_ps_annot_en_adj, axis=0)
    diffsq = np.square(all_ps_annot_en_adj-mean_all_ps_annot_en_adj)
    annot_en_adj_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return annot_en_adj, annot_en_adj_se, all_ps_annot_en_adj

def get_gcorsq_enrichment(gcorsq, gcorsq_se, ps_gcorsq, use_jk_adj):
    """
    Estimate enrichment of genetic correlation
    """
    
    # get dimensions
    nannot = gcorsq.shape[0]
    nblock = ps_gcorsq.shape[0]

    # get all-sample estimate
    all_gcorsq_en = gcorsq/gcorsq[0]
   
    # get variance of bottom
    gcorsq_var = np.square(gcorsq_se)
    var_btm = gcorsq_var[0]

    # get covariance of top and bottom
    nblock, nannot = ps_gcorsq.shape
    cov_top_btm_tmp = np.zeros((nblock, nannot))
    for i in range(nannot):
        cov_top_btm_tmp[:,i] = (ps_gcorsq[:,i]-np.mean(ps_gcorsq[:,i])) * \
            (ps_gcorsq[:,0]-np.mean(ps_gcorsq[:,0]))
    cov_top_btm = np.mean(cov_top_btm_tmp, axis=0)
    cov_top_btm *= (np.float64(nblock)-1.0)
  
    # adjust for bias
    if use_jk_adj == False:
        all_gcorsq_en_adj = (all_gcorsq_en + cov_top_btm/gcorsq[0]) / \
            (1.0 + var_btm/gcorsq[0])
    else:
        # no adjustment if using jackknife bias correction
        # pseudo estimate likely not reliable
        all_gcorsq_en_adj = all_gcorsq_en.copy()

    # get jackknife pseudo estimates
    all_ps_gcorsq_en = np.zeros((nblock, nannot), dtype=np.float64)
    factor = np.float64(nblock)/(np.float64(nblock)-1.0)
    for i in range(nblock):
        annot_est = ps_gcorsq[i,:]/ps_gcorsq[i,0]
       
        # adjust for bias if using analytical bias correction
        if use_jk_adj == False:
            cov_top_btm_tmp = factor*cov_top_btm
            var_btm_tmp = factor*var_btm
            anont_est = (annot_est + cov_top_btm_tmp/ps_gcorsq[i,0]) / \
                (1.0 + var_btm_tmp/ps_gcorsq[i,0])
            all_ps_gcorsq_en[i,:] = annot_est
        else:
            all_ps_gcorsq_en[i,:] = annot_est

    # get se
    nblock = np.float64(nblock)
    mean_all_ps_gcorsq_en = np.mean(all_ps_gcorsq_en, axis=0)
    diffsq = np.square(all_ps_gcorsq_en-mean_all_ps_gcorsq_en)
    all_gcorsq_en_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return all_gcorsq_en, all_gcorsq_en_se

# get genetic correlation estimate
def get_gcor(hsq1, ps_hsq1, hsq2, ps_hsq2, gcov, ps_gcov, bound):
    
    # get estimate
    nannot = hsq1.shape[0]
    annot_est = np.zeros(nannot, dtype=np.float32)
    idx = (hsq1>0) & (hsq2>0)
    annot_est[idx] = gcov[idx] / np.sqrt(hsq1[idx]*hsq2[idx])
    if bound == True:
        annot_est[annot_est>1.0] = 1.0
        annot_est[annot_est<-1.0] = -1.0

    # get standard error
    nblock = ps_hsq1.shape[0]
    all_ps_annot_est = np.zeros((nblock, nannot))
    for j in range(nblock):
        ps_annot_est = np.zeros(nannot, dtype=np.float32)
        idx = (ps_hsq1[j,:]>0) & (ps_hsq2[j,:]>0)
        top = ps_gcov[j,idx]
        btm = np.sqrt(ps_hsq1[j,idx]*ps_hsq2[j,idx])
        ps_annot_est[idx]= top / btm
        if bound == True:
            ps_annot_est[ps_annot_est>1.0] = 1.0
            ps_annot_est[ps_annot_est<-1.0] = -1.0
        all_ps_annot_est[j,:] = ps_annot_est

    # get se
    nblock = np.float64(nblock)
    mean_all_ps_annot_est = np.mean(all_ps_annot_est, axis=0)
    diffsq = np.square(all_ps_annot_est-mean_all_ps_annot_est)
    annot_est_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return annot_est, annot_est_se

def shrink_sum(prior, prior_var, ps_prior, m_prior,
    post, post_var, ps_post, m_post, factor):
   
    nparam = post.shape[0]
    nblock = ps_post.shape[0]

    shrunk = (factor*post/m_post + (1-factor)*prior/m_prior) * m_post
    ps_shrunk = np.zeros((nblock, nparam))
    for i in range(nblock):
        ps_shrunk[i,:] = \
          (factor*ps_post[i,:]/m_post+(1-factor)*ps_prior[i]/m_prior)*m_post

    diffsq = np.square(ps_shrunk-shrunk)
    shrunk_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return shrunk, shrunk_se, ps_shrunk

def get_shrink_factor(est, est_var, annot_nsnp, shrinkage):

    prior_var = est_var[0]/annot_nsnp[0]
    post_var = est_var/annot_nsnp
    factor = 1.0 / (1.0 + shrinkage*post_var/prior_var)

    return factor

def get_gcorsq(hsq1, ps_hsq1, hsq2, ps_hsq2, gcov, ps_gcov, annot_nsnp,
    bound, shrinkage, use_jk_adj):

    """
    Estimate squared genetic correlation
    """
  
    # get dimensions
    nblock = ps_hsq1.shape[0]
    nannot = hsq1.shape[0]
   
    # apply shrinkage on heritability and genetic covariance
    hsq1_var = (nblock-1)*np.mean(np.square(ps_hsq1-hsq1), axis=0)
    hsq2_var = (nblock-1)*np.mean(np.square(ps_hsq2-hsq2), axis=0)
    gcov_var = (nblock-1)*np.mean(np.square(ps_gcov-gcov), axis=0)
    factor1 = get_shrink_factor(hsq1, hsq1_var, annot_nsnp, shrinkage)
    factor2 = get_shrink_factor(hsq2, hsq2_var, annot_nsnp, shrinkage)
    factorx = get_shrink_factor(gcov, gcov_var, annot_nsnp, shrinkage)
    factor = np.fmin(np.fmin(factor1, factor2), factorx)
    hsq1_shrink, _, ps_hsq1_shrink = shrink_sum(hsq1[0], hsq1_var[0],
        ps_hsq1[:,0], annot_nsnp[0], hsq1, hsq1_var, ps_hsq1,
        annot_nsnp, factor)
    hsq2_shrink, _, ps_hsq2_shrink = shrink_sum(hsq2[0], hsq2_var[0],
        ps_hsq2[:,0], annot_nsnp[0], hsq2, hsq2_var, ps_hsq2,
        annot_nsnp, factor)
    gcov_shrink, _, ps_gcov_shrink = shrink_sum(gcov[0], gcov_var[0],
        ps_gcov[:,0], annot_nsnp[0], gcov, gcov_var, ps_gcov,
        annot_nsnp, factor)

    # get estimates
    gcov_var_shrink = get_jkcov(gcov_shrink, ps_gcov_shrink,
        gcov_shrink, ps_gcov_shrink)
    hsq_jkcov_shrink = get_jkcov(hsq1_shrink, ps_hsq1_shrink,
        hsq2_shrink, ps_hsq2_shrink)
    top = np.square(gcov_shrink) - gcov_var_shrink
    btm = hsq1_shrink*hsq2_shrink - hsq_jkcov_shrink
    annot_est = top/btm
   
    # get pseudo estimates
    factor = np.float64(nblock)/(np.float64(nblock)-1.0)
    ps_top = np.square(ps_gcov_shrink) - factor*gcov_var_shrink
    ps_btm = ps_hsq1_shrink*ps_hsq2_shrink - factor*hsq_jkcov_shrink

    # get variance of bottom covariance of top and bottom
    var_btm = get_jkcov(btm, ps_btm, btm, ps_btm)
    cov_top_btm = get_jkcov(top, ps_top, btm, ps_btm)
    
    # adjust for bias
    if use_jk_adj == False:
        annot_est_adj = (annot_est + cov_top_btm/btm) / (1.0 + var_btm/btm)

    # apply bound
    if bound == True:
        annot_est[annot_est>1.0] = 1.0
        annot_est[annot_est<0.0] = 0.0

    # get pseudo values
    all_ps_annot_est = np.zeros((nblock, nannot), dtype=np.float64)
    factor = np.float64(nblock)/(np.float64(nblock)-1.0)
    for j in range(nblock):
        
        # get estimate
        topj = ps_top[j,:]
        btmj = ps_btm[j,:]
        ps_annot_est = topj/btmj

        # apply bound
        if bound == True:
            ps_annot_est[ps_annot_est>1.0] = 1.0
            ps_annot_est[ps_annot_est<0.0] = 0.0
 
        # adjust for bias
        if use_jk_adj == False:
            var_btm_tmp = factor*var_btm
            cov_top_btm_tmp = factor*cov_top_btm
            ps_annot_est = (ps_annot_est + cov_top_btm_tmp/btmj) / \
                (1.0 + var_btm_tmp/btmj)

        # store pseudo estimate
        all_ps_annot_est[j,:] = ps_annot_est

    # adjust for bias using jackknife
    if use_jk_adj == True:
        jknife_mean = np.mean(all_ps_annot_est, axis=0)
        bias = (nblock - 1) * (jknife_mean - annot_est)
        annot_est_adj = annot_est - bias
        factor = np.float64(nblock)/(np.float64(nblock)-1.0)
        all_ps_annot_est = all_ps_annot_est - factor*bias

    # get se
    mean_all_ps_annot_est = np.mean(all_ps_annot_est, axis=0)
    diffsq = np.square(all_ps_annot_est-mean_all_ps_annot_est)
    annot_est_adj_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return annot_est_adj.astype(np.float64), \
           annot_est_adj_se.astype(np.float64), \
           all_ps_annot_est.astype(np.float64)

def get_gcovsq_diff(hsq1, ps_hsq1, hsq2, ps_hsq2, gcov, gcov_se,
    ps_gcov, gcorsq, ps_gcorsq):

    """
    Estimate squared genetic correlation
    """
  
    # get cov of hsq
    hsq_jkcov = get_jkcov(hsq1, ps_hsq1, hsq2, ps_hsq2)

    # get dimensions
    nblock = ps_hsq1.shape[0]
    nannot = hsq1.shape[0]

    # get estimates
    top = np.square(gcov) - np.square(gcov_se)
    btm = gcorsq[0]*(hsq1*hsq2 - hsq_jkcov)
    annot_est = top - btm

    # base should be exactly 0
    annot_est[0] = 0.0

    # get pseudo values
    all_ps_annot_est = np.zeros((nblock, nannot), dtype=np.float64)
    for j in range(nblock):
        
        gcovj = ps_gcov[j,:]
        gcorsqj = ps_gcorsq[j,:]
        hsq1j = ps_hsq1[j,:]
        hsq2j = ps_hsq2[j,:]

        factor = np.float64(nblock)/(np.float64(nblock)-1.0)
        top = np.square(gcovj) - factor*np.square(gcov_se)
        btm = gcorsqj[0]*(hsq1j*hsq2j - factor*hsq_jkcov)

        ps_annot_est = top - btm
        ps_annot_est[0] = 0.0
        all_ps_annot_est[j,:] = ps_annot_est

    # get se
    nblock = np.float64(nblock)
    mean_all_ps_annot_est = np.mean(all_ps_annot_est, axis=0)
    diffsq = np.square(all_ps_annot_est-mean_all_ps_annot_est)
    annot_est_se = np.sqrt((nblock-1)*np.mean(diffsq, axis=0))

    return annot_est, annot_est_se, all_ps_annot_est

def get_weights(ldscore1, ldscore2):

    ldscore1_pos = np.fmax(ldscore1, 1.0)
    ldscore2_pos = np.fmax(ldscore2, 1.0)
    ldscorex_pos = np.sqrt(ldscore1_pos*ldscore2_pos)

    w1 = 1.0/ldscore1_pos
    w2 = 1.0/ldscore2_pos
    wx = 1.0/ldscorex_pos

    return w1, w2, wx

def get_pred(coef, x, n1, n2, intercept):
    
    nprod = np.sqrt(n1*n2)
    score_sum = np.sum(x[:,:-1], axis=1)
    pred = coef*score_sum*nprod + intercept

    return pred.astype(np.float32)

def update_weights(w1, w2, wx, pred1, pred2, predx):
   
    var1 = 2.0*np.square(pred1)
    var2 = 2.0*np.square(pred2)
    var12 = np.sqrt(var1*var2)/2.0 + np.square(predx)

    new_w1 = w1 / var1
    new_w2 = w2 / var2
    new_wx = wx / var12

    return new_w1, new_w2, new_wx

def load_frqfile(frqfile_fnm, start_chrom, stop_chrom):
    
    all_frq = []
    for i in range(start_chrom, stop_chrom+1):
        all_frq.append(pd.read_table('{}{}.frq'.format(frqfile_fnm,i),
            delim_whitespace=True))
    
    all_frq = pd.concat(all_frq, axis=0, ignore_index=True)
    sigma = np.sqrt(2.0*all_frq['MAF']*(1.0-all_frq['MAF']))
    all_frq['SIGMA'] = sigma.astype(np.float64)

    return all_frq

def estimate_gcor(sumstats_fnm, score_fnm, weight_fnm, annot_fnm,
    frqfile_fnm, out_fnm, fit_intercept, use_intercept, save_ps_coef,
    use_chrom, use_robust_regression, nblk_jk, min_maf, bound,
    apply_shrinkage, use_jk_adj):
   
    """
    Estimate coefficients, enrichments, and squared genetic correlations
    """

    # load summary stats file and perform filtering
    sumstats1 = load_sumstats(sumstats_fnm[0])
    sumstats2 = load_sumstats(sumstats_fnm[1])

    # filter summary stats
    sumstats1 = filter_sumstats(sumstats1)
    sumstats2 = filter_sumstats(sumstats2)

    # load scores
    start_chrom, stop_chrom = use_chrom
    score_snps, annot_names, ldscore1 = load_score(score_fnm, 'pop1',
        start_chrom, stop_chrom)
    _, _, ldscore2 = load_score(score_fnm, 'pop2', start_chrom, stop_chrom)
    _, _, ldscorex = load_score(score_fnm, 'te', start_chrom, stop_chrom)
    
    logging.info('Loaded LD scores for {} SNPs {} annotations'\
        .format(score_snps.shape[0], annot_names.shape[0]))

    # load annotation matrix
    annot_snps, _, annot_mat = load_annot(annot_fnm, start_chrom, stop_chrom)
    
    logging.info('Loaded annotations: {} SNPs, {} annotations'\
        .format(annot_mat.shape[0], annot_mat.shape[1]))

    # load regression snps
    regsnps, _, ldscore1_reg = load_score(weight_fnm, 'pop1',
        start_chrom, stop_chrom)
    
    _,_,ldscore2_reg = load_score(weight_fnm, 'pop2', start_chrom, stop_chrom)

    logging.info('Loaded LD scores for regression SNPs')

    # load frqfile
    freq1 = load_frqfile(frqfile_fnm[0], start_chrom, stop_chrom)
    freq2 = load_frqfile(frqfile_fnm[1], start_chrom, stop_chrom)
    
    logging.info('Loaded allele frequency files')

    # intersect snps
    snp1_idx,snp2_idx,score_idx,reg_idx,freq_idx = intersect_snp_regression(
        sumstats1['SNP'], sumstats2['SNP'], score_snps['SNP'],
        regsnps['SNP'], freq1['SNP'])
    
    sumstats1 = sumstats1.loc[snp1_idx,:]
    sumstats2 = sumstats2.loc[snp2_idx,:]
    
    ldscore1 = ldscore1[score_idx,:]
    ldscore2 = ldscore2[score_idx,:]
    ldscorex = ldscorex[score_idx,:]
    
    ldscore1_reg = ldscore1_reg[reg_idx,0]
    ldscore2_reg = ldscore2_reg[reg_idx,0]
    
    sigma1 = freq1.loc[freq_idx,'SIGMA'].values.astype(np.float32)
    sigma2 = freq2.loc[freq_idx,'SIGMA'].values.astype(np.float32)

    nregsnp = ldscore1.shape[0]

    logging.info('After intersection, {} SNPs are left for regression'\
        .format(sumstats1.shape[0]))

    # get regression weights
    weight1,weight2,weightx = get_weights(ldscore1_reg, ldscore2_reg)
    
    # estimate initial coefficient
    tau1,int1 = get_coef_raw(ldscore1, sumstats1, sumstats1, 1.0)
    tau2,int2 = get_coef_raw(ldscore2, sumstats2, sumstats2, 1.0)
    theta,intx = get_coef_raw(ldscorex, sumstats1, sumstats2, 0.0)

    # predict chi-square and prodcut
    n1,n2 = sumstats1['N'].values,sumstats2['N'].values
    pred1 = get_pred(tau1, ldscore1, n1, n1, int1)
    pred2 = get_pred(tau2, ldscore2, n2, n2, int2)
    predx = get_pred(theta, ldscorex, n1, n2, intx)

    # update weight
    weight1_, weight2_, weightx_ = update_weights(weight1, weight2, weightx,
        pred1, pred2, predx)
    
    # subtraction when intercept is constrained
    subtract1 = np.ones(nregsnp, dtype=np.float32)*use_intercept[0]
    subtract2 = np.ones(nregsnp, dtype=np.float32)*use_intercept[1]
    subtractx = np.ones(nregsnp, dtype=np.float32)*use_intercept[2]

    # get regression coefficients
    tau1, tau1_se, ps_tau1 = get_coef(ldscore1, sumstats1, sumstats1,
       weight1_, nblk_jk, fit_intercept[0], use_robust_regression, subtract1)
    tau2, tau2_se, ps_tau2 = get_coef(ldscore2, sumstats2, sumstats2,
        weight2_, nblk_jk, fit_intercept[1], use_robust_regression, subtract2)
    theta,theta_se,ps_theta = get_coef(ldscorex, sumstats1, sumstats2,
        weightx_, nblk_jk, fit_intercept[2], use_robust_regression, subtractx)

    logging.info('Obtained coefficients')

    # save pseudo coefficient if user asks
    if save_ps_coef == True:
        np.savetxt('{}.pseudo_tau1.gz'.format(out_fnm), ps_tau1)
        np.savetxt('{}.pseudo_tau2.gz'.format(out_fnm), ps_tau2)
        np.savetxt('{}.pseudo_theta.gz'.format(out_fnm), ps_theta)
        
        logging.info('Saved pseudo coefficients')

    # get snps for estimation
    annot_idx = intersect_snp_estimation(annot_snps['SNP'],freq1,freq2,min_maf)
    annot_mat = annot_mat[annot_idx,:]

    # get hsq and gcov for each annotation
    hsq1, hsq1_se, ps_hsq1 = get_sum(annot_mat, tau1, ps_tau1)
    hsq2, hsq2_se, ps_hsq2 = get_sum(annot_mat, tau2, ps_tau2)
    gcov, gcov_se, ps_gcov = get_sum(annot_mat, theta, ps_theta)
   
    # get enrichment estimate
    annot_nsnp, annot_std = get_annot_sumstat(annot_mat)
    hsq1_en,hsq1_en_se,ps_hsq1_en = get_enrichment(hsq1, ps_hsq1, annot_nsnp)
    hsq2_en,hsq2_en_se,ps_hsq2_en = get_enrichment(hsq2, ps_hsq2, annot_nsnp)
    gcov_en,gcov_en_se,ps_gcov_en = get_enrichment(gcov, ps_gcov, annot_nsnp)

    # get genetic correlation
    gcor,gcor_se = get_gcor(hsq1, ps_hsq1, hsq2, ps_hsq2, gcov, ps_gcov, bound)

    # get sqaured genetic correlation estimates
    gcorsq, gcorsq_se, ps_gcorsq = get_gcorsq(hsq1, ps_hsq1, hsq2, ps_hsq2,
        gcov, ps_gcov, annot_nsnp, bound, apply_shrinkage, use_jk_adj)

    # estimate enrichment of squared genetic correlation
    gcorsq_en, gcorsq_en_se = get_gcorsq_enrichment(gcorsq,
        gcorsq_se, ps_gcorsq, use_jk_adj)
    tstat = np.fabs((gcorsq_en-1.0) / (gcorsq_en_se+1e-16))
    gcorsq_en_pval = (1.0-scipy.stats.t.cdf(tstat, nblk_jk-1))*2.0

    # difference in gcov squared
    gcovsq_diff, gcovsq_diff_se, ps_gcovsq_diff = get_gcovsq_diff(hsq1,
        ps_hsq1, hsq2, ps_hsq2, gcov, gcov_se, ps_gcov, gcorsq, ps_gcorsq)
    tstat = np.fabs(gcovsq_diff / (gcovsq_diff_se+1e-16))
    gcovsq_diff_pval = (1.0-scipy.stats.t.cdf(tstat, nblk_jk-1))*2.0

    # collect results in a data frame
    out_df = pd.DataFrame()
    out_df['ANNOT'] = annot_names['ANNOT']
    out_df['NSNP'] = annot_nsnp
    out_df['STD'] = annot_std
    out_df['TAU1'] = tau1[:-1]; out_df['TAU1_SE'] = tau1_se[:-1]
    out_df['TAU2'] = tau2[:-1]; out_df['TAU2_SE'] = tau2_se[:-1]
    out_df['THETA'] = theta[:-1]; out_df['THETA_SE'] = theta_se[:-1]
    out_df['HSQ1'] = hsq1; out_df['HSQ1_SE'] = hsq1_se
    out_df['HSQ2'] = hsq2; out_df['HSQ2_SE'] = hsq2_se
    out_df['GCOV'] = gcov; out_df['GCOV_SE'] = gcov_se
    out_df['GCOR'] = gcor; out_df['GCOR_SE'] = gcor_se
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
    out_df.to_csv(out_fnm, sep='\t', index=False)

    # print out other related information
    logging.info('hsq1: {} {}'.format(hsq1[0], hsq1_se[0]))
    if fit_intercept[0] == 'yes':
        logging.info('intercept hsq1: {} {}'.format(tau1[-1], tau1_se[-1]))
    else:
        logging.info('intercept hsq1: {} {}'.format(use_intercept[0], 0.0))
    logging.info('hsq2: {} {}'.format(hsq2[0], hsq2_se[0]))
    if fit_intercept[1] == 'yes':
        logging.info('intercept hsq2: {} {}'.format(tau2[-1], tau2_se[-1]))
    else:
        logging.info('intercept hsq2: {} {}'.format(use_intercept[1], 0.0))
    logging.info('gcov: {} {}'.format(gcov[0], gcov_se[0]))
    if fit_intercept[2] == 'yes':
        logging.info('intercept gcov: {} {}'.format(theta[-1], theta_se[-1]))
    else:
        logging.info('intercept gcov: {} {}'.format(use_intercept[2], 0.0))
    logging.info('gcor: {} {}'.format(gcor[0], gcor_se[0]))
    logging.info('gcorsq: {} {}'.format(gcorsq[0], gcorsq_se[0]))
