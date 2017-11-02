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
