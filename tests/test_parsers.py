"""
Test the filter functions
"""

import os
from os import path
import unittest
import logging
import pprint
import sys
import json
import csv

from operator import itemgetter
from itertools import ifilter
from collections import namedtuple, defaultdict
from msings.utils import walker
from msings import parsers, filters


from __init__ import TestBase
import __init__ as config

log = logging.getLogger(__name__)

testMSIfile = path.join(config.datadir, 'MSI')

class TestParsers(TestBase):
    """Test each of the parsers are returning the correct fieldnames, 
    prefixes and variant_key list"""
    def setUp(self):
        self.outdir = self.mkoutdir()

    def tearDown(self):
        pass

    def testMSIParser(self):
        specimens = defaultdict(dict)
        prefixes = []
        variant_keys = []
        control_info=open(path.join(testMSIfile, 'testMSIcontrol'),'rU')
        files = walker(testMSIfile)        
        files = filter(filters.msi_file_finder,files)
        analysis_type='parsers.parse_msi'
        multiplier=2.0
        threshold=[0.2, 0.9]
        chosen_parser='{}(files, control_info, specimens, prefixes, variant_keys, multiplier,threshold)'.format(analysis_type)    
        specimens, prefixes, fieldnames, variant_keys=eval(chosen_parser)  
        self.assertListEqual(sorted(prefixes),sorted(['0228T_CON_OPXv4_INT', '5437_E05_OPXv4_NA12878_MA0013', '6037_E05_OPXv4_NA12878_HA0201']))
        self.assertListEqual(sorted(fieldnames), sorted(['0228T_CON_OPXv4_INT', '5437_E05_OPXv4_NA12878_MA0013', '6037_E05_OPXv4_NA12878_HA0201', 'Position']))
        self.assertListEqual(variant_keys, ['Position'])
        
