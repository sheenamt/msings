"""
Create the summary file of all MSI_Analysis files in folder, requires output of count_msi_samples for each sample. No recalculation happening in this script.

Usage:

msi create_summary $SAVEPATH -o $OUTFILE

"""
import logging
import csv
import sys
import argparse
import collections
import os
from itertools import ifilter
from msings import filters
from msings.utils import walker
import pandas as pd

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('path',
                        help='Path to analysis files')
    parser.add_argument('-o','--outfile', type = argparse.FileType('w'),
                        default = sys.stdout,
                        help='Name of the output file')
    

def action(args):
    top_level_df=pd.DataFrame()
    msi_files = ifilter(filters.msi_analysis,walker(args.path))
    msi_files=sorted(msi_files)    

    for pth in msi_files:
        df_pfxinfo=pd.read_csv(os.path.join(pth.dir, pth.fname), delimiter='\t')
        if top_level_df.empty:
            top_level_df=df_pfxinfo
        else:
            top_level_df=pd.merge(top_level_df, df_pfxinfo, how='left')

    top_level_df.to_csv(args.outfile, na_rep=' ', index=False, header=True, sep='\t')
