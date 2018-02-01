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
from msings import parsers


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
        analysis_type='parsers.parse_msi'
        multiplier=2.0
        threshold=[0.2, 0.9]
        chosen_parser='{}(files, control_info, specimens, prefixes, variant_keys, multiplier,threshold)'.format(analysis_type)    
        specimens, prefixes, variant_keys=eval(chosen_parser)  

        variant_keys=[variant_keys]
        fieldnames = variant_keys + list(prefixes) 
        self.assertListEqual(sorted(prefixes),sorted(['0228T', '5437_NA12878', '6037_NA12878']))
        self.assertListEqual(sorted(fieldnames), sorted(['0228T', '5437_NA12878', '6037_NA12878', 'Position']))
        self.assertListEqual(variant_keys, ['Position'])
        
