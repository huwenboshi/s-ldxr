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

    # get analysis code
    analysis, argmap = check_command_line(args)

    # run the corresponding analysis
    if analysis == 'estimate_score':
        
        estimate_score(argmap['annot'], argmap['score'], argmap['bfile'],
            argmap['ld-wind-cm'], argmap['print-snps'], argmap['out'])

    elif analysis == 'estimate_gcor':
        
        estimate_gcor(argmap['gcor'], argmap['ref-ld-chr'], argmap['w-ld-chr'],
            argmap['annot'], argmap['frqfile'], argmap['out'],
            argmap['add-intercept'], argmap['save-pseudo-coef'],
            argmap['use-chrom'], argmap['use-robust-regression'],
            argmap['n-blocks'], argmap['min-maf'], argmap['bound'],
            argmap['apply-shrinkage'])

    else:
        logging.error('Invalid command line option.')

    # end the log
    end_log()


# check command line
def check_command_line(args):

    # parse command line
    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    # for estimating score
    if(argmap['score'] != None and argmap['annot'] != None and
       argmap['ld-wind-cm'] != None and argmap['bfile'] != None and
       argmap['print-snps'] != None and argmap['out'] != None):
        
        if argmap['score'] != 'standardized' and argmap['score'] != 'allelic':
            return ('invalid', argmap)

        return ('estimate_score', argmap)

    # for performing regression and estimation
    if(argmap['gcor'] != None and argmap['out'] != None and
       argmap['ref-ld-chr'] != None and
       argmap['frqfile'] != None and argmap['annot'] and
       argmap['save-pseudo-coef'] != None and
       argmap['add-intercept'] != None and argmap['w-ld-chr'] != None and
       argmap['use-chrom'] != None and argmap['n-blocks'] != None and
       argmap['use-robust-regression'] !=None and argmap['min-maf'] != None
       and argmap['bound'] != None and argmap['apply-shrinkage'] != None):
        
        return ('estimate_gcor', argmap)

    # command line option not recognized
    return ('invalid', argmap)

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

# get command line input
def get_command_line():
    
    # create the help document
    parser = argparse.ArgumentParser(description='Partition trans-ethnic '
        'genetic correlation by functional categories')

    ########## estimate score ##########

    parser.add_argument('--score', dest='score', default='standardized',
        required=False, type=str, help='Estimate score')

    parser.add_argument('--bfile', dest='bfile', type=str, required=False,
        default=None, nargs=2, help='Reference panel file in PLINK format')

    parser.add_argument('--annot', dest='annot', type=str, required=False,
        nargs='+', help='Annotation file')

    parser.add_argument('--ld-wind-cm', dest='ld-wind-cm', type=float,
        required=False, help='Window size in cm', default=1.0)

    parser.add_argument('--print-snps', dest='print-snps', type=str,
        required=False, help='A list of snps to output the scores')

    ########## estimate genetic correlation ##########

    parser.add_argument('--gcor', dest='gcor', type=str, nargs=2,
        required=False, help='GWAS summary statitistics')

    parser.add_argument('--ref-ld-chr', dest='ref-ld-chr', type=str,
        required=False, nargs='+', help='Prefix for the score files')

    parser.add_argument('--w-ld-chr', dest='w-ld-chr', type=str,
        required=False, nargs='+', help='Prefix for the weight files')

    parser.add_argument('--add-intercept', dest='add-intercept', type=str,
        required=False, default=['yes','yes','no'], nargs=3,
        help='Specific whether to use intercept in regression')

    parser.add_argument('--use-chrom', dest='use-chrom', type=int,
        required=False, default=[1, 22], nargs=2,
        help='Specific which chromosomes to use')

    parser.add_argument('--n-blocks', dest='n-blocks', type=int,
        required=False, default=200, help='Number of jackknife blocks')

    parser.add_argument('--use-robust-regression',dest='use-robust-regression',
        default=False, required=False, action='store_true',
        help='Use robust regression instead of OLS')

    parser.add_argument('--frqfile', dest='frqfile', type=str, nargs=2,
        required=False, help='Prefix of the frequency files')

    parser.add_argument('--min-maf', dest='min-maf', type=float, default=0.05,
        required=False, help='Minimum MAF for performing the estimation')

    parser.add_argument('--save-pseudo-coef', dest='save-pseudo-coef',
        default=False, required=False, action='store_true',
        help='Save jack knife pseudo coefficient')

    parser.add_argument('--bound', dest='bound',
        default=False, required=False, action='store_true',
        help='Bound the genetic correlation estimate')

    parser.add_argument('--apply-shrinkage', dest='apply-shrinkage',
        default=0.5, required=False, type=float,
        help='Apply shrinkage on the estimates')

    ########## output  ##########

    # specify output file
    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

# execute the main function
if(__name__ == '__main__'):
    main()
