import logging, gzip, sys
import numpy as np
import pandas as pd
from .annot import *
from pysnptools.snpreader import Bed

maf_thres = 0.001
eps = 1e-8

def load_legend(bfile_prefix):
    """
    Load legend of the plink file
    """
    
    legend = pd.read_table(bfile_prefix+'.bim',
        delim_whitespace=True, header=None)
    
    legend.columns = ['CHR', 'SNP', 'CM', 'BP', 'A2', 'A1']
    
    return legend

def load_refpanel(bfile, load_idx):
    """
    Load the genotype data
    """

    # read data
    genotype = bfile[:,load_idx].read(dtype=np.float32).val
    
    # impute with average for missing snps
    nanidx = np.where(np.isnan(genotype))
    mean_geno = np.nanmean(genotype, axis=0, dtype=np.float32)
    genotype[nanidx] = mean_geno[nanidx[1]]
    nindv = genotype.shape[0]

    # get minor allele frequency
    maf = np.sum(genotype, axis=0, dtype=np.float32) / (2.0*np.float32(nindv))
    maf[maf>0.5] = 1.0-maf[maf>0.5]
    mac = maf*2*np.float32(nindv)

    # get sigmasq
    sigmasq = genotype.var(axis=0, dtype=np.float32) + eps

    # center genotype matrix
    genotype -= genotype.mean(axis=0, dtype=np.float32)
    
    return genotype, nindv, maf, mac, sigmasq

def get_ld_allelic(snpdata, start, stop):
    """
    Get LD matrix of a region of genotype data
    """

    nindv = np.float32(snpdata.shape[0])
    ld = np.dot(snpdata[:,start:stop].T, snpdata)/(nindv-1.0)
  
    return ld

def get_ld_standardized(snpdata, start, stop, sigmasq):
    """
    Get LD matrix of a region of genotype data
    """

    nindv = np.float32(snpdata.shape[0])
    ld = np.dot(snpdata[:,start:stop].T, snpdata)/(nindv-1.0)
  
    ld = ld / np.sqrt(sigmasq)
    ld = (ld.T / np.sqrt(sigmasq[start:stop])).T
    np.fill_diagonal(ld[:,start:stop], 1)

    return ld

def get_ldscore_mat_standardized(ld1, n1, maf1, ld2, n2, maf2, start, stop):
    """
    Estimate unbiased squared LD matirx and matrix products
    """

    # compute square and product the ld matrices
    ld1_sq = np.square(ld1)
    ld2_sq = np.square(ld2)
    ld1_ld2 = np.multiply(ld1, ld2)
    
    # adjust squared ld matrix (no adjustment for product)
    ld1_sq -= (1.0-ld1_sq) / (np.float32(n1)-2.0)
    ld2_sq -= (1.0-ld2_sq) / (np.float32(n2)-2.0)

    # zero out entries involving snps with maf less than maf_thres
    ld1_sq[maf1[start:stop]<=maf_thres,:] = 0.0
    ld2_sq[maf2[start:stop]<=maf_thres,:] = 0.0
    ld1_ld2[(maf1[start:stop]<=maf_thres)|(maf2[start:stop]<=maf_thres),:]=0.0
    
    ld1_sq[:,maf1<=maf_thres] = 0.0
    ld2_sq[:,maf2<=maf_thres] = 0.0
    ld1_ld2[:,(maf1<=maf_thres)|(maf2<maf_thres)] = 0.0

    # fill diagonal with 1
    np.fill_diagonal(ld1_sq[:,start:stop], 1.0)
    np.fill_diagonal(ld2_sq[:,start:stop], 1.0)
    np.fill_diagonal(ld1_ld2[:,start:stop], 1.0)

    return ld1_sq, ld2_sq, ld1_ld2

def get_ldscore_mat_allelic(ld1, n1, maf1, sigmasq1, ld2, n2, maf2,
    sigmasq2, start, stop):
    """
    Estimate unbiased squared LD matirx and matrix products
    """

    # compute square and product the ld matrices
    ld1_sq = np.square(ld1, dtype=np.float32)
    ld2_sq = np.square(ld2, dtype=np.float32)
    ld1_ld2 = np.multiply(ld1, ld2)
    
    # adjust squared ld matrix (no adjustment needed for product)
    n1,n2 = np.float32(n1),np.float32(n2)
    sigmasq1_sub = sigmasq1[start:stop]
    sigmasq2_sub = sigmasq2[start:stop]
    ld1_sq -= np.outer(sigmasq1_sub, sigmasq1) / np.float32(n1-1.0)
    ld2_sq -= np.outer(sigmasq2_sub, sigmasq2) / np.float32(n2-1.0)
    ld1_sq *= n1/(n1-1.0); ld2_sq *= n2/(n2-1.0)

    # zero out entries involving snps with maf less than maf_thres
    ld1_sq[maf1[start:stop]<=maf_thres,:] = 0.0
    ld2_sq[maf2[start:stop]<=maf_thres,:] = 0.0
    ld1_ld2[(maf1[start:stop]<=maf_thres)|(maf2[start:stop]<=maf_thres),:]=0.0
    
    ld1_sq[:,maf1<=maf_thres] = 0.0
    ld2_sq[:,maf2<=maf_thres] = 0.0
    ld1_ld2[:,(maf1<=maf_thres)|(maf2<maf_thres)] = 0.0

    # fill diagonal with sigmasq
    np.fill_diagonal(ld1_sq[:,start:stop],
        (n1-1.0)/(n1+1.0)*np.square(sigmasq1_sub))
    np.fill_diagonal(ld2_sq[:,start:stop],
        (n2-1.0)/(n2+1.0)*np.square(sigmasq2_sub))
    np.fill_diagonal(ld1_ld2[:,start:stop], sigmasq1_sub*sigmasq2_sub)

    # divide by sigmasq
    sigma12_sub = np.sqrt(sigmasq1_sub*sigmasq2_sub)
    ld1_sq = ld1_sq / sigmasq1_sub[:,np.newaxis]
    ld2_sq = ld2_sq / sigmasq2_sub[:,np.newaxis]
    ld1_ld2 = ld1_ld2 / sigma12_sub[:,np.newaxis]

    return ld1_sq, ld2_sq, ld1_ld2

# calculate the scores
def calc_scores(score_type, legend, bfile1, bfile2, annot, win, annot_names):
  
    # get data infor
    tot_nsnp = legend.shape[0]
    tot_nannot = annot_names.shape[0]

    # initialize scores
    ldsc1 = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)
    ldsc2 = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)
    ldscx = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)

    # extract snp cm information
    all_snp_cm = legend['CM'].values

    # process nsnp_per_iter snps at a time, iterate until all snps traversed
    start_idx = 0
    nsnp_per_iter = 2000
    while True:
        
        # from start (inclusive) to stop (exclusive)
        stop_idx = start_idx + nsnp_per_iter
        if stop_idx > tot_nsnp:
            stop_idx = tot_nsnp

        # get legend for snps to be processed
        snp_cm = all_snp_cm[start_idx:stop_idx]

        # find the range of snps to load
        start_cm = all_snp_cm[start_idx] - win
        stop_cm = all_snp_cm[stop_idx-1] + win

        load_idx = np.where((all_snp_cm>=start_cm) &
                            (all_snp_cm<=stop_cm))[0]
        
        # find index of start_idx and (stop_idx-1) within load_idx
        start = np.where(load_idx==start_idx)[0][0]
        stop = np.where(load_idx==(stop_idx-1))[0][0]+1
    
        # load centered genotype data
        geno1, nref1, maf1, mac1, sigmasq1 = load_refpanel(bfile1, load_idx)
        geno2, nref2, maf2, mac2, sigmasq2 = load_refpanel(bfile2, load_idx)

        # get squared / product of ld matrices for standardized genotypes
        if score_type == 'standardized':
            ld1 = get_ld_standardized(geno1, start, stop, sigmasq1)
            ld2 = get_ld_standardized(geno2, start, stop, sigmasq2)
            ld1_sq,ld2_sq,ld1_ld2 = get_ldscore_mat_standardized(ld1, nref1,
                maf1, ld2, nref2, maf2, start, stop)
        
        # get squared / product of ld matrices for centered genotypes
        if score_type == 'allelic':
            ld1 = get_ld_allelic(geno1, start, stop)
            ld2 = get_ld_allelic(geno2, start, stop)
            
            ld1_sq,ld2_sq,ld1_ld2 = get_ldscore_mat_allelic(ld1, nref1, maf1,
                sigmasq1, ld2, nref2, maf2, sigmasq2, start, stop)
    
        # compute ld score matrix
        nsnp = snp_cm.shape[0]
        annot_load = annot[load_idx,:]
        load_cm = all_snp_cm[load_idx]
        ldsc1_mat = ldsc1[start_idx:stop_idx,:]
        ldsc2_mat = ldsc2[start_idx:stop_idx,:]
        ldscx_mat = ldscx[start_idx:stop_idx,:]
        for i in range(nsnp):
            cm = snp_cm[i]
            tmp_idx = np.where((load_cm>=(cm-win)) & (load_cm<=(cm+win)))[0]
            tmp_mat = annot_load[tmp_idx,:]
            ldsc1_mat[i,:] = np.dot(ld1_sq[i,tmp_idx], tmp_mat)
            ldsc2_mat[i,:] = np.dot(ld2_sq[i,tmp_idx], tmp_mat)
            ldscx_mat[i,:] = np.dot(ld1_ld2[i,tmp_idx], tmp_mat)

        # check stop condition
        if stop_idx == tot_nsnp:
            break
        
        # update start idx
        start_idx = stop_idx

    return ldsc1, ldsc2, ldscx

def write_score(all_annot, score_mat, legend, printsnps, out_fnm):
    """
    Write score to file
    """
    out_f = gzip.open(out_fnm, 'w')

    # write header
    header_fields = ['CHR', 'SNP', 'BP'] + all_annot['ANNOT'].values.tolist()
    header = '\t'.join(header_fields)
    header += '\n'
    out_f.write(header.encode())

    all_chr = legend['CHR'].values
    all_snp = legend['SNP'].values
    all_bp = legend['BP'].values
    nsnp = legend.shape[0]
    for i in range(nsnp):
        if all_snp[i] not in printsnps:
            continue
        line = '{}\t{}\t{}\t'.format(all_chr[i], all_snp[i], all_bp[i])
        line += '\t'.join([str(c) for c in score_mat[i,:].tolist()])
        line += '\n'
        out_f.write(line.encode())

    out_f.close()

def estimate_score(annot_fnm, score_type, bfile_prefix, wind,
    printsnps, out_fnm):
    """
    Estimate annotation specific socres for an annotation file
    """

    # load annotation matrix
    all_snp, all_annot, annot_mat = load_annot_chrom(annot_fnm)

    # regression snps
    snps = pd.read_table(printsnps, header=None)

    # load legend data
    legend1 = load_legend(bfile_prefix[0])
    legend2 = load_legend(bfile_prefix[1])

    # estimate scores
    bfile1 = Bed(bfile_prefix[0], count_A1=False)
    bfile2 = Bed(bfile_prefix[1], count_A1=False)
    ldsc1, ldsc2, ldscx = calc_scores(score_type, legend1, bfile1, bfile2,
        annot_mat, wind, all_annot)

    # write scores to file
    snp_set = set(snps[0].tolist())
    write_score(all_annot,ldsc1,legend1,snp_set,'{}_pop1.gz'.format(out_fnm))
    write_score(all_annot,ldsc2,legend2,snp_set,'{}_pop2.gz'.format(out_fnm))
    write_score(all_annot,ldscx,legend1,snp_set,'{}_te.gz'.format(out_fnm))


def load_score(prefix_list, suffix, start_chrom, end_chrom,
    annot_start_idx=3, snp_idx=1):
    """
    Load scores for all chromosomes
    """

    # get snps
    all_snp = []
    for i in range(start_chrom, end_chrom+1):
        prefix = prefix_list[0]
        filename = '{}{}_{}.gz'.format(prefix, i, suffix)
        tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
            na_filter=False, memory_map=True, usecols=['SNP'])
        all_snp.append(tbl)
    all_snp = pd.concat(all_snp, axis=0, ignore_index=True)
    tot_nsnp = all_snp.shape[0]

    # get annotations
    all_annot = []
    annot_list = []
    for i in range(len(prefix_list)):
        prefix = prefix_list[i]
        filename = '{}{}_{}.gz'.format(prefix, start_chrom, suffix)
        with gzip.open(filename, 'r') as f:
            line = f.readline().strip().decode("utf-8")
            tmp = line.split()[annot_start_idx:]
            all_annot += tmp
            annot_list.append(tmp)
    tot_nannot = len(all_annot)
    
    # load the scores into memory, extra column for storing intercept
    all_score = np.zeros((tot_nsnp, tot_nannot+1), dtype=np.float32)
    all_score[:,-1] = np.float32(1.0)
    idx_c = 0
    for k in range(len(prefix_list)):
        prefix = prefix_list[k]
        idx_r = 0
        for i in range(start_chrom, end_chrom+1):
            filename = '{}{}_{}.gz'.format(prefix, i, suffix)
            tmp = annot_list[k]
            dt_load = dict(zip(tmp, [np.float32]*len(tmp)))
            tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
                na_filter=False, memory_map=True, usecols=tmp, dtype=dt_load)
            nrow, ncol = tbl.shape
            all_score[idx_r:idx_r+nrow,idx_c:idx_c+ncol] = tbl[tmp].values
            idx_r += nrow
        idx_c += ncol

    # convert all annot to data frame
    all_annot = pd.DataFrame({'ANNOT': all_annot})

    return all_snp, all_annot, all_score
