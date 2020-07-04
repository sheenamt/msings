"""
Create MSI Baseline file from multiple MSI negative samples

Usage:

msi create_baseline  /path/to/sample/files -o output_file

"""
import os
import csv
import sys
import argparse
import re
import natsort

from collections import defaultdict, namedtuple
from itertools import groupby
from operator import itemgetter
from numpy import std, array, average
from msings.utils import walker
from msings import filters

def build_parser(parser):

    parser.add_argument('path', 
                        help='Path to analysis files')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')
    

def action(args):
    control = defaultdict(list)
    # apply a series of filters to files
    files = filter(filters.msi_file_finder, walker(args.path))
    #sort the files so that the output in the workbook is sorted
    files = natsort.natsorted(files)

    for pth in files:
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            info = sorted(reader, key=itemgetter('Position'))
            for k, g in groupby(info, key=itemgetter('Position')):
                for row in g:
                    if int(row['Average_Depth']) >= 30:
                        control[row['Position']].append(int(row['Number_of_Peaks']))
    header = ['Position', 'Standard_Deviation', 'Average_Number_Peaks', 'Count']
    writer = csv.writer(args.outfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
    writer.writerow(header)
    for k, v in sorted(control.items()):
        count = len(v)
        a = array(v)
        std_d = a.std()
        ave = average(a)
        row = [k, std_d, ave, count]
        if count >=3:
            writer.writerow(row)

