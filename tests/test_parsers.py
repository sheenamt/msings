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
import pandas as pd

from pandas.util.testing import assert_frame_equal # <-- for testing dataframes
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

    def testPasrseThreshold(self):
        """ Test the the max and min are correct, given 1 or 2 inputs, error on 2+ inputs"""
        self.assertEqual((0.2, 0.2), parsers.parse_thresholds([0.2]))
        self.assertEqual((0.1, 0.2), parsers.parse_thresholds([0.2, 0.1]))
        self.assertEqual((0.2, 0.9), parsers.parse_thresholds([0.2, 0.9]))
        self.assertRaises(ValueError, parsers.parse_thresholds, [0.2, 0.9, 0.3])

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
        

    def testTumorMutationBurden(self):
        specimens = pd.DataFrame()
        d = [{'Position':'tumor_mutation_burden','0228T':'6/3006'}]
        expectedDF = pd.DataFrame(data=d)
        prefixes = []
        files = walker(testMSIfile)
        specimens = parsers.parse_total_mutation_burden(specimens, prefixes, files)
        #ensure columns are the same
        self.assertListEqual(sorted(expectedDF.columns), sorted(specimens.columns))
        #Ensure 'position' was added'
        self.assertItemsEqual(sorted(expectedDF["Position"]), sorted(specimens["Position"]))
        #Ensure calculation is correct
        self.assertItemsEqual(sorted(expectedDF["0228T"]), sorted(specimens["0228T"]))

    def testFilter(self):
        """Test the OPX/BRO filter"""
        variant_to_include=['frameshift insertion', 'frameshift deletion', 'frameshift block substitution', 'stopgain', 
                             'stoploss', 'nonframeshift insertion', 'nonframeshift deletion', 
                             'nonframeshift block substitution', 'nonsynonymous SNV', 'synonymous SNV', 'exonic', 'splicing']
        uw_freq_P={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        uw_freq_F={'UW_Freq': '0.5','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        exac_P={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        exac_F={'UW_Freq': '0.0005','EXAC': '0.5',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        var_reads_P={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic, intronic'}
        var_reads_F={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '4','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        thouG_P={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        thouG_F={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '0.5', 'Variant_Type': 'frameshift deletion,exonic'}
        var_type_P={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'frameshift deletion,exonic'}
        var_type_F={'UW_Freq': '0.0005','EXAC': '-1',  'Var_Reads': '8','1000g_ALL': '-1', 'Variant_Type': 'intronic'}

        self.assertTrue(parsers.opx_bro_filter(variant_to_include, uw_freq_P))
        self.assertFalse(parsers.opx_bro_filter(variant_to_include, uw_freq_F))
        self.assertTrue(parsers.opx_bro_filter(variant_to_include, exac_P))
        self.assertFalse(parsers.opx_bro_filter(variant_to_include, exac_F))
        self.assertTrue(parsers.opx_bro_filter(variant_to_include, var_reads_P))
        self.assertFalse(parsers.opx_bro_filter(variant_to_include, var_reads_F))
        self.assertTrue(parsers.opx_bro_filter(variant_to_include, thouG_P))
        self.assertFalse(parsers.opx_bro_filter(variant_to_include, thouG_F))
        self.assertTrue(parsers.opx_bro_filter(variant_to_include, var_type_P))
        self.assertFalse(parsers.opx_bro_filter(variant_to_include, var_type_F))





