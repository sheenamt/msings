"""
Test the utils functions
"""

import os
from os import path
import unittest
import logging
import pprint
import sys
import json

from msings.utils import munge_pfx, munge_path, munge_date

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

class TestUtils(TestBase):

    def setUp(self):
        self.run=munge_path('/home/genetics/data/130321_HA00211_OncoPlex64')
        self.outdir = self.mkoutdir()

    def testMungePFX(self):
        real_info={'control': 'NA12878', 
                   'machine-run': 'HA0201', 
                   'library-version': 'OPXv4', 
                   'well': 'E05', 
                   'run': '60', 
                   'sample_id': '6037', 
                   'pfx': '6037_E05_OPXv4_NA12878_HA0201',
                   'assay':'OncoPlex',
                   'mini-pfx': '6037_NA12878'}
        
        test_info=munge_pfx('6037_E05_OPXv4_NA12878_HA0201')
        self.assertDictEqual(real_info, test_info)

