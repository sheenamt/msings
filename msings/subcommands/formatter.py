"""
Formats MSI bed file for input into varscan's readcount function.

Usage:

msi formatter /path/to/bed/file -o /path/to/outfile

"""
import os
import csv
import sys
import argparse
import re
import natsort
from collections import defaultdict

def build_parser(parser):
    parser.add_argument('bedfile', 
                        type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to tab delimited bed file of format chr start end')
    parser.add_argument('outfile', 
                        type=argparse.FileType('w'),
                        help='Name of the output file')

def coords(row, chromosome=0, start=1, stop=2):
    """Return chromosome, start, end - assuming that column indices are
    the same in bedfile and giantfile.

    """

    return row[chromosome], int(row[start]), int(row[stop])

def msi_interval_creator(ranges):
    """Change the range information into varscan format
    chr start - T 
    """
    data=[]
    #make sure the chr is sorted
    for chr in natsort.natsorted(ranges.keys()):
        #make sure the start position is sorted
        for row in sorted(ranges[chr]):
            info=(chr, row, '-', 'T')
            data.append(info)
    return data

def action(args):
    bedfile=csv.reader(args.bedfile,delimiter='\t')
    output=args.outfile
        
    # prepare a dictionary of chromosome: set(positions)
    # includes all positions between start-stop
    ranges = defaultdict(set)
    msi_calls=defaultdict()
    writer=csv.writer(output,delimiter='\t')
    for row in bedfile:
        chr, beg, end = coords(row)
        ranges[chr].update(range(beg, end + 1))
    writer.writerows(msi_interval_creator(ranges))







