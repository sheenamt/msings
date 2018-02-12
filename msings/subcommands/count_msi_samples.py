"""
Create msi output file with MSI location counts for multiple samples

Usage:

msi count_msi_samples  /path/to/baseline/file /path/to/sample/files -m 2.5 -t 0.1 0.7 -o output_file

"""
import os
import csv
import sys
import argparse
import re
import natsort

from collections import defaultdict, namedtuple
from itertools import groupby, ifilter
from operator import itemgetter
from numpy import std, array, average, sum

import pandas as pd

from msings.parsers import parse_msi, parse_total_mutation_burden
from msings.utils import walker

def build_parser(parser):
    parser.add_argument('control_file', 
                        type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to control file')
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-m','--multiplier',
                        default=2.0,
                        type=float,
                        help='STD cutoff multiplier, default of 2.0')
    parser.add_argument('-t','--msi_threshold',
                        nargs='+',
                        help='MSI score threshold or range of two thresholds, default of 0.2')
    parser.add_argument('-b', '--tumor_burden', 
                        action='store_true',
                        help='Calcuatel tumor burden. FOR INTERNAL UW USE ONLY, WILL FAIL OTHERWISE')
    parser.add_argument('-o', '--outfile', 
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')

def action(args):
        
    specimens = defaultdict(dict)
    prefixes = []
    variant_keys =[]
    files = list(walker(args.path))
    analysis_type='parse_msi'

    multiplier=args.multiplier
    if args.msi_threshold:
        threshold=args.msi_threshold        
    else:
        threshold=[0.2,0.2]
    control_file = args.control_file
    chosen_parser='{}(files, control_file, specimens, prefixes, variant_keys, multiplier,threshold)'.format(analysis_type)
    specimens, prefixes, variant_keys =eval(chosen_parser)

    #Column headers for output 
    fieldnames = [variant_keys] + list(prefixes) 

    #list of fields to print first in output
    msi_fields =['unstable_loci', 'passing_loci','msings_score','msi_status'] 

    #only run tumor_burden at UW
    if args.tumor_burden:
        specimens=parse_total_mutation_burden(specimens, prefixes, files)
        msi_fields.append('tumor_mutation_burden')

    writer = csv.writer(args.outfile, delimiter = '\t')
    writer.writerow(fieldnames)
    
    #next print the msi status info, then remove from dataframe before printing each loci detail
    for info in msi_fields:
       parsed_info = specimens['Position']==info
       to_print=specimens[parsed_info].to_csv(args.outfile, sep='\t',header=False, index=False, columns=fieldnames, na_rep = ' ')
       #Drop that row from specimens
       specimens = specimens.ix[~(specimens['Position']==info)]

    specimens.to_csv(args.outfile, na_rep=' ', index=False, columns=fieldnames, header=False, sep='\t')

