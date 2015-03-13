import os
import csv
import sys
import re
import argparse

from collections import defaultdict, namedtuple
from itertools import groupby, ifilter
from operator import itemgetter
from numpy import std, array, average

parser = argparse.ArgumentParser()
parser.add_argument('control_file', type=argparse.FileType('rU'),
                    default=sys.stdin,
                    help='Path to control file')
parser.add_argument('path',
                    help='Path to analysis files')
parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                    default=sys.stdout,
                    help='Name of the output file')
args = parser.parse_args()

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
    prefixes = []
    control_info = csv.DictReader(args.control_file, delimiter="\t")
    sample_info = defaultdict(list)
    # apply a series of filters to files
    files = ifilter(msi_file_finder, walker(args.path))
    #sort the files so that the output in the workbook is sorted
    files = sorted(files)

    for row in control_info:
        for pth in files:
            pfx = pth.fname.split('.')[0]
            if pfx not in prefixes:
                prefixes.append(pfx)
            with open(os.path.join(args.path, pth.fname)) as fname:
                reader = csv.DictReader(fname, delimiter='\t')
                reader = sorted(reader, key=itemgetter('Position'))
                for key, group in groupby(reader, key=itemgetter('Position')):
                    #for each position in control (row[0])
                    if row['Position'] == key:
                        for info in group:
                            if int(info['Avg_read_depth']) >= 30:
                                value = float(row['Ave']) + (2 * float(row['Std']))
                                if int(info['Number_Peaks']) >= value:
                                    new_info = 1
                                else:
                                    new_info = 0
                            else:
                                new_info = 'NA'
                            sample_info[row['Position']].append(new_info)

    annotation_headers = ['Position']
    annotation_headers.extend(prefixes)
    writer = csv.writer(args.outfile, delimiter='\t')
    writer.writerow(annotation_headers)
    for k, v in sorted(sample_info.items()):
#        print k, v
        line = [k]
        line.extend(v)
        writer.writerow(line)



if __name__=='__main__':
    action(args)
