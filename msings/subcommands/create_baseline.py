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
from itertools import groupby, ifilter
from operator import itemgetter
from numpy import std, array, average, sum

from msings.parsers import parse_msi
from msings.utils import walker

def build_parser(parser):

    parser.add_argument('path', 
                        help='Path to analysis files')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Name of the output file')
    

def msi_file_finder(pth):
    """
    Return True if pth represents an analysis file.
    """
    return bool(re.search(r'.msi.txt', pth.fname))
    

Path = namedtuple('Path', ['dir','fname'])

def walker(dir):
    """Recursively traverse direcory `dir`. For each tuple containing
    (path, dirs, files) return a named tuple with attributes (dir,
    fname) for each file name in `files`.
    """
    for (pth, dirs, files) in os.walk(dir):
        for fname in files:
            yield Path(pth, fname)

def action(args):
    control = defaultdict(list)
    # apply a series of filters to files
    files = ifilter(msi_file_finder, walker(args.path))
    #sort the files so that the output in the workbook is sorted
    files = sorted(files)

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
    for k, v in natsort.natsorted(control.items()):
        count = len(v)
        a = array(v)
        std_d = a.std()
        ave = average(a)
        row = [k, std_d, ave, count]
        if count >=3:
            writer.writerow(row)

