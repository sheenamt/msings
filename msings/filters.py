"""
Filters for analysis files. In each function, pth is a namedtuple
with attributes (dir, fname).
"""

import re

def msi_file_finder(pth):
    """
    Return True if pth represents a msi file file.
    """
    return bool(pth.fname.endswith('.msi.txt'))

def snp_analysis(pth):
    """
    True only for pfx.Analysis.{txt,csv}
    """
    return bool(pth.fname.endswith('SNP_Analysis.txt'))

def msi_analysis(pth):
    """
    True only for pfx.Analysis.{txt,csv}
    """
    return bool(pth.fname.endswith('MSI_Analysis.txt'))
