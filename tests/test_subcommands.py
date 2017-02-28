"""
Test the subcommand scripts
"""

import os
from os import path
import unittest
import logging
import csv
import sys
import json
import copy

from numpy import std, average
from operator import itemgetter
from itertools import groupby

from msings.subcommands import analyzer
from msings.subcommands import count_msi_samples
from msings.subcommands import create_baseline
from msings.subcommands import formatter


from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)

msi_testfiles = path.join(config.datadir, 'MSI')

control ='5437_E05_OPXv4_NA12878_MA0013'

MSI_LOCI={'1': {1: '1:1-5', 2: '1:1-5', 3: '1:1-5', 4: '1:1-5', 5: '1:1-5', 7: '1:7-11', 8: '1:7-11', 9: '1:7-11', 10: '1:7-11', 11: '1:7-11'}, 
          '7': {1: '7:1-5', 2: '7:1-5', 3: '7:1-5', 4: '7:1-5', 5: '7:1-5'}}

OUTPUT_RAW={'1:1-5': {'total_depth': 0, 'wildtype_depth': 0, 'Name': 'NAME1', 'mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}, 
            '7:1-5': {'total_depth': 0, 'wildtype_depth': 0, 'Name': 'NAME3', 'mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}, 
            '1:7-11': {'total_depth': 0, 'wildtype_depth': 0, 'Name': 'NAME2', 'mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}}

MSI_SITE_DATA={'1:1-5': {'site_depth': 100, 'total_depth': 500, 'wildtype_depth': 500, 'indels': {},
                         'mutant_depth': 0,'mutant_tally': 0, 'total_sites': 5, 'Name': 'NAME1'}, 
               '7:1-5': {'site_depth': 50, 'total_depth': 250, 'wildtype_depth': 136,
                         'indels': {1: {'site_depth': 150, 'mutant_tally': 3, 'allele_fraction': 0.7333333333333333, 'mutant_depth': 110}, 
                                    -1: {'site_depth': 50, 'mutant_tally': 1, 'allele_fraction': 0.08, 'mutant_depth': 4}}, 
                         'mutant_depth': 114, 'mutant_tally': 4, 'total_sites': 5, 'Name': 'NAME3'}, 
               '1:7-11': {'site_depth': 100, 'total_depth': 500, 'wildtype_depth': 210,
                          'indels': {1: {'site_depth': 300, 'mutant_tally': 3, 'allele_fraction': 0.4, 'mutant_depth': 120}, 
                                     -3: {'site_depth': 300, 'mutant_tally': 3, 'allele_fraction': 0.35, 'mutant_depth': 105}}, 
                          'mutant_depth': 225, 'mutant_tally': 6, 'total_sites': 5, 'Name': 'NAME2'}}

OUTPUT={'1:1-5': {'IndelLength:AlleleFraction:Reads': '0:1.0:100', 'Standard_Deviation': 0, 'Average_Depth': 100, 'Number_of_Peaks': 1, 'Name': 'NAME1'}, 
        '7:1-5': {'IndelLength:AlleleFraction:Reads': '-1:0.109090909091:4 0:0.741818181818:27 1:1.0:110', 'Standard_Deviation': '0.493303', 'Average_Depth': 50, 'Number_of_Peaks': 3, 'Name': 'NAME3'}, 
        '1:7-11': {'IndelLength:AlleleFraction:Reads': '-3:0.875:105 -2:0:0 -1:0:0 0:0.0:27 1:1.0:120', 'Standard_Deviation': '1.904576', 'Average_Depth': 100, 'Number_of_Peaks': 3, 'Name': 'NAME2'}}

class TestFormatter(TestBase):
    """
    Test the msi formatter script
    """
    def testCoords(self):
        
        row =['1', '45795895', '45795905', 'NAME1']
        output =('1', 45795895, 45795905)
        self.assertEqual(formatter.coords(row), output)

    def testMSIIntervalCreator(self):
        ranges={'1': set([45795904, 45795905, 45795895, 45795896])}
        data=[('1', 45795895, '-', 'T'), ('1', 45795896, '-', 'T'),
              ('1', 45795904, '-', 'T'), ('1', 45795905, '-', 'T')]

        self.assertEqual(formatter.msi_interval_creator(ranges), data)

class TestAnalyzer(TestBase):
    """
    Test the msi analyzer subcommands
    """

    def testParseMSIBedfile(self):
        """Test that the MSI Bed file is parsed correctly
        """
        msi_sites, output_info={}, {}
        for row in csv.DictReader(open(path.join(msi_testfiles, 'test.msi.bed')), delimiter='\t', fieldnames=['chrom','start','end','name']):
            msi_sites, output_info = analyzer.parse_msi_bedfile(row, msi_sites, output_info)
        self.assertDictEqual(msi_sites, MSI_LOCI)
        self.assertDictEqual(output_info, OUTPUT_RAW)
        
    def testCalcMSIDist(self):
        """Test MSI site distribution calculation"""
        self.maxDiff = None
        output_info = copy.deepcopy(OUTPUT_RAW)
        sample_msi=csv.DictReader(open(path.join(msi_testfiles, 'test.msi_output')), delimiter='\t', restkey='Misc')

        for row in sample_msi:
            loci_position = MSI_LOCI[row['chrom']][int(row['position'])]
            output_info[loci_position].update(analyzer.calc_msi_dist(row, output_info[loci_position]))

        self.assertDictEqual(output_info, MSI_SITE_DATA)
        
    def testCalcSummaryStats(self):
        """Test MSI summary calculations
        """
        self.maxDiff=None
        local_msi_site = copy.deepcopy(MSI_SITE_DATA)
        output_local={}
        cutoff=float(0.05)
        output_local.update(analyzer.calc_summary_stats(local_msi_site, cutoff))
        self.assertDictEqual(output_local, OUTPUT)               

    def testHighestPeak(self):
        """Test that the highest peak is returned
        """
        msi_sites1=copy.deepcopy(MSI_SITE_DATA['7:1-5'])
        msi_sites2=copy.deepcopy(MSI_SITE_DATA['1:7-11'])
        wt_1 = float(msi_sites1['wildtype_depth'])/msi_sites1['total_depth']
        wt_2 = float(msi_sites2['wildtype_depth'])/msi_sites2['total_depth']
        highest_peak1 = analyzer.calc_highest_peak(msi_sites1['indels'], wt_1)
        highest_peak2 = analyzer.calc_highest_peak(msi_sites2['indels'], wt_2)
        self.assertEqual(0.7333333333333333, highest_peak1)
        self.assertEqual(0.42, highest_peak2)
        
    def testCalcNumberPeaks(self):
        """Test that the number of peaks and the peak annotation
        is being calculated/parsed correctly. 
        """
        msi_sites1=copy.deepcopy(MSI_SITE_DATA['1:7-11'])
        wt_1 = float(msi_sites1['wildtype_depth'])/msi_sites1['total_depth']
        cutoff=[0.05]
        peaks = []
        sites={0:'0:0:0'}
        peaks, sites=analyzer.calc_number_peaks(msi_sites1['indels'], sites, 0.42, cutoff)
        output_peaks=1
        output_site_info= {0: '0:0:0', 1: '1:0.952380952381:120', -3: '-3:0.833333333333:105'}
        self.assertEqual(peaks, output_peaks)
        self.assertEqual(sites, output_site_info)

    def testCalcWildType(self):
        """Test the Wildtype calculations"""
        msi_sites=copy.deepcopy(MSI_SITE_DATA['1:7-11'])
        sites = {}

        wt_1 = float(msi_sites['wildtype_depth'])/msi_sites['total_depth']
        wt_ave = float(msi_sites['wildtype_depth'])/msi_sites['total_sites']
        sites=analyzer.calc_wildtype(msi_sites['indels'].keys(), wt_ave, wt_1, 0.42)
        wt_output={0: '0:1.0:42.0', -1: '-1:0:0', -3: '-3:0:0', -2: '-2:0:0', 1: '1:0:0'}

        self.assertDictEqual(sites, wt_output)

    def testCalcSTDPeaks(self):
        """Test the standard deviation calculations"""
        peaks=['0:0.863414634146:354', '1:0.0402598525993:17', '-1:0.0855382887727:34', '-2:0.0132135895294:5']
        stdev=analyzer.calc_std_peaks(peaks)
        self.assertEqual(stdev, '0.410894')

