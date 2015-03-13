"""
Compare sample MSI loci count to the specified control file

Usage:

msings msi_sample_vs_control /path/to/control/msi/file /path/to/sample/msi/file -o output_file

"""

import os
import csv
import sys
import argparse
import collections

from itertools import groupby, ifilter
from operator import itemgetter


def build_parser(parser):
    parser.add_argument('control_file', type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to control file')
    parser.add_argument('sample_files', nargs='+',
                        help='Path to sample files')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')


def action(args):
    specimens = defaultdict(dict)
    prefixes = []
    variant_keys =[]
    files = walker(args.path)  
    analysis_type='parse_msi'

    if args.msithreshold:
        threshold=args.msithreshold
    else:
        threshold=[0.2,]
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




