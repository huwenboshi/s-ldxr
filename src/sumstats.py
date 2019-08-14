import numpy as np
import pandas as pd
import os, sys, gzip, logging

# required columns for the summary stats file
required = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N']

# magic bits to discern file type
magic_dict = {
    "\x1f\x8b\x08": "gz",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip"
    }
max_len = max(len(x) for x in magic_dict)

# define equivalent alleles
equiv = dict()
equiv["AC"] = set(["TG", "AC", "TC", "AG"])
equiv["AG"] = set(["TC", "AG", "TG", "AC"])
equiv["CA"] = set(["GT", "CA", "GA", "CT"])
equiv["CT"] = set(["GA", "CT", "GT", "CA"])
equiv["TC"] = set(["AG", "TC", "AC", "TG"])
equiv["TG"] = set(["AC", "TG", "AG", "TC"])
equiv["GA"] = set(["CT", "GA", "CA", "GT"])
equiv["GT"] = set(["CA", "GT", "CT", "GA"])

# define reversed alleles
reverse = dict()
reverse["AC"] = set(["GT", "CA", "CT", "GA"])
reverse["AG"] = set(["CT", "GA", "GT", "CA"])
reverse["CA"] = set(["TG", "AC", "AG", "TC"])
reverse["CT"] = set(["AG", "TC", "TG", "AC"])
reverse["TC"] = set(["GA", "CT", "CA", "GT"])
reverse["TG"] = set(["CA", "GT", "GA", "CT"])
reverse["GA"] = set(["TC", "AG", "AC", "TG"])
reverse["GT"] = set(["AC", "TG", "TC", "AG"])

# define strand ambiguous alleles
ambiguous = set(["AT", "CG", "TA", "GC"])

def load_sumstats(filename):
    """
    Load GWAS summary association.
    Perform initial filtering, including removing SNPs without rs ID and
    SNPs with allele length greater than 1.

    Args:
        filename: File name of the summary association data.

    Returns:
        A data frame containing the GWAS summary association data.

    Raises:
        KeyError: Raises an exception.
    """
        
    # load summary stats
    sumstats_cols = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N']
    sumstats = pd.read_table(filename, delim_whitespace=True, engine='c',
        na_filter=False, memory_map=True, usecols=sumstats_cols,
        dtype={'Z': np.float32, 'N': np.float32}) 

    # log the results
    logging.info('Loaded {} SNPs from the GWAS summary data file'\
        .format(sumstats.shape[0]))

    return sumstats

def filter_sumstats(sumstats):
    """
    Filter out SNPs with ambiguous rs ID, BP, alleles.
    """
        
    # remove duplicates based on both rs ID and position
    sumstats = sumstats.drop_duplicates('SNP', keep=False)
    #sumstats = sumstats.drop_duplicates('BP', keep=False)
    sumstats = sumstats.reset_index(drop=True)

    # remove snps with extremely large chis square
    chisq_max = max(0.001*sumstats['N'].max(), 80)
    chisq = np.square(sumstats['Z'])
    sumstats = sumstats.drop(np.where(chisq>chisq_max)[0])
    sumstats = sumstats.reset_index(drop=True)

    logging.info('{} SNPs left after filtering'\
        .format(sumstats.shape[0]))

    return sumstats
