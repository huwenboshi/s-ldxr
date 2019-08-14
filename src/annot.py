import logging, gzip, sys, mmap
import numpy as np
import pandas as pd

def load_annot_chrom(filename_list, annot_start_idx=4, snp_idx=2):
    """
    Load annotation of a single chromosome
    """

    # get snps snps
    all_snp = pd.read_table(filename_list[0], delim_whitespace=True,
        engine='c', na_filter=False, memory_map=True, usecols=['SNP'])
    tot_nsnp = all_snp.shape[0]

    # count number of annotations
    all_annot = []
    annot_list = []
    for i in xrange(len(filename_list)):
        filename = filename_list[i]
        with gzip.open(filename, 'r') as f:
            line = f.readline().strip()
            cols = line.split()[annot_start_idx:]
            all_annot += cols
            annot_list.append(cols)
    tot_nannot = len(all_annot)
    
    # load the annot into memory
    all_annot_mat = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)
    idx_c = 0
    for i in xrange(len(filename_list)):
        filename = filename_list[i]
        tmp = annot_list[i]
        dt_load = dict(zip(tmp, [np.float32]*len(tmp)))
        tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
            na_filter=False, memory_map=True, usecols=tmp, dtype=dt_load)
        ncol = tbl.shape[1]
        all_annot_mat[:,idx_c:idx_c+ncol] = tbl[tmp].values
        idx_c += ncol

    # convert annotation to data frame
    all_annot = pd.DataFrame({'ANNOT': all_annot})

    return all_snp, all_annot, all_annot_mat

def load_annot(prefix_list, start_chrom, end_chrom,
    annot_start_idx=4, snp_idx=2):
    """
    Load annotation of all chromosomes
    """

    # load all snps snps
    all_snp = []
    for i in xrange(start_chrom, end_chrom+1):
        prefix = prefix_list[0]
        filename = '{}{}.annot.gz'.format(prefix, i)
        tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
            na_filter=False, memory_map=True, usecols=['SNP'])
        all_snp.append(tbl)
    all_snp = pd.concat(all_snp, axis=0, ignore_index=True)
    tot_nsnp = all_snp.shape[0]

    # get annotations
    all_annot = []
    annot_list = []
    for i in xrange(len(prefix_list)):
        prefix = prefix_list[i]
        filename = '{}{}.annot.gz'.format(prefix, start_chrom)
        with gzip.open(filename) as f:
            line = f.readline().strip()
            tmp = line.split()[annot_start_idx:]
            all_annot += tmp
            annot_list.append(tmp)
    tot_nannot = len(all_annot)
    
    # load the annot into memory
    all_annot_mat = np.zeros((tot_nsnp, tot_nannot), dtype=np.float32)
    idx_c = 0
    for k in xrange(len(prefix_list)):
        prefix = prefix_list[k]
        idx_r = 0
        for i in xrange(start_chrom, end_chrom+1):
            filename = '{}{}.annot.gz'.format(prefix, i)
            tmp = annot_list[k]
            dt_load = dict(zip(tmp, [np.float32]*len(tmp)))
            tbl = pd.read_table(filename, delim_whitespace=True, engine='c',
                na_filter=False, memory_map=True, usecols=tmp, dtype=dt_load)
            nrow, ncol = tbl.shape
            all_annot_mat[idx_r:idx_r+nrow,idx_c:idx_c+ncol] = tbl[tmp].values
            idx_r += nrow
        idx_c += ncol
    
    # convert annot to data frame
    all_annot = pd.DataFrame({'ANNOT': all_annot})

    return all_snp, all_annot, all_annot_mat

def get_annot_sumstat(annot_mat):
    """
    Get summary statistics of the annotation matrix
    """
    annot_nsnp = np.sum(annot_mat, axis=0)
    annot_std = np.std(annot_mat, axis=0)
    #annot_corr = np.corrcoef(annot_mat.T)

    return annot_nsnp, annot_std #, annot_corr
