"""
Create msi output file with MSI location counts, can be used to create output for one or many samples

Usage:

msi count_msi_samples  /path/to/baseline/file /path/to/sample/file -m 2.5 -t 0.1 0.7 -o output_file

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
from msings import filters

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
    parser.add_argument('-p', '--pipeline_manifest', type=argparse.FileType('rU'),
                        help='Path to pipeline manifest, used for ordering output in UW pipeline specifically')    
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
    files_sorted=[]
    if args.pipeline_manifest:
        sort_order = [x['barcode_id'] for x in csv.DictReader(args.pipeline_manifest)]
        for sample in sort_order:
            #Grab the file for each sample, in specified sort order
            pfx_file = [s for s in files if sample in s.fname]
            if pfx_file:
                files_sorted.append(pfx_file[0])
    else:
        files_sorted=files

    analysis_type='parse_msi'
    multiplier=args.multiplier
    if args.msi_threshold:
        threshold=args.msi_threshold        
    else:
        threshold=[0.2,0.2]
    control_file = args.control_file

    chosen_parser='{}(files, control_file, specimens, prefixes, variant_keys, multiplier,threshold)'.format(analysis_type)
    df_specimens, prefixes, variant_keys =eval(chosen_parser)

    #Column headers for output 
    fieldnames = [variant_keys] + list(prefixes) 

    #list of fields to print first in output
    msi_fields =['unstable_loci', 'covered_loci','msings_score','msi_status'] 

    #only run tumor_burden at UW
    if args.tumor_burden:
        df_specimens=parse_total_mutation_burden(df_specimens, prefixes, files)
        msi_fields.append('tumor_mutation_burden')

    writer = csv.writer(args.outfile,delimiter = '\t')
    fieldnames.strip()
    writer.writerow(fieldnames)
    
    #next print the msi status info, then remove from dataframe before printing each loci detail
    for info in msi_fields:
       parsed_info = df_specimens['Position']==info
       to_print=df_specimens[parsed_info].to_csv(args.outfile, sep='\t',header=False, index=False, columns=fieldnames, na_rep = ' ')
       #Drop that row from specimens
       df_specimens = df_specimens.ix[~(df_specimens['Position']==info)]

    df_specimens.to_csv(args.outfile, na_rep=' ', index=False, columns=fieldnames, header=False, sep='\t')


