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

from collections import defaultdict, namedtuple
from itertools import groupby, ifilter
from operator import itemgetter
from numpy import std, array, average, sum

from msings.parsers import parse_msi
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
    parser.add_argument('-o', '--outfile', 
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')

def action(args):
        
    specimens = defaultdict(dict)
    prefixes = []
    variant_keys =[]
    files = walker(args.path)  
    analysis_type='parse_msi'

    multiplier=args.multiplier
    if args.msi_threshold:
        threshold=args.msi_threshold        
    else:
        threshold=[0.2,0.2]
    control_file = args.control_file
    chosen_parser='{}(files, control_file, specimens, prefixes, variant_keys, multiplier,threshold)'.format(analysis_type)
    specimens, prefixes, fieldnames, variant_keys=eval(chosen_parser)

    writer = csv.DictWriter(args.outfile, fieldnames = fieldnames,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    #for position in sorted(specimens.keys(), reverse=True):
    for position in sorted(specimens.keys(), reverse=True):
        d = {k:v for k,v in zip(variant_keys,position)}  
        d.update({pfx:specimens[position].get(pfx) for pfx in prefixes})  
        writer.writerow(d)


