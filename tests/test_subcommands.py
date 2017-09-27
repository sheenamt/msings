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

from numpy import std, average,ceil
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
          '7': {1: '7:1-5', 2: '7:1-5', 3: '7:1-5', 4: '7:1-5', 5: '7:1-5', 7: '7:7-11', 8: '7:7-11', 9: '7:7-11', 10: '7:7-11', 11: '7:7-11'},
          '8': {1: '8:1-5', 2: '8:1-5', 3: '8:1-5', 4: '8:1-5', 5: '8:1-5'}}

OUTPUT_RAW={'1:1-5': {'total_depth': 0, 'Name': 'WT-ONLY', 'total_mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}, 
            '7:1-5': {'total_depth': 0, 'Name': 'MUT-BIG>AVE', 'total_mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}, 
            '1:7-11': {'total_depth': 0, 'Name': 'WT-BIGGEST', 'total_mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}},
            '7:7-11': {'total_depth': 0,  'Name': 'NO-COV', 'total_mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}},
            '8:1-5': {'total_depth': 0, 'Name': 'MUT-BIG<AVE', 'total_mutant_depth': 0, 'mutant_tally': 0, 'total_sites': 0, 'indels': {}}}

MSI_SITE_DATA={'1:1-5': {'site_depth': 100, 'total_depth': 500, 'Name': 'WT-ONLY', 'mutant_tally': 0, 'total_mutant_depth': 0, 'total_sites': 5, 'indels': {}}, 
               '1:7-11': {'site_depth': 100, 'total_depth': 500, 'Name': 'WT-BIGGEST', 'mutant_tally': 6, 'total_mutant_depth': 30, 'total_sites': 5, 
                          'indels': {1: {'site_depth': 300, 'mutant_tally': 3, 'allele_fraction': 0, 'mutant_depth': 14}, -3: {'site_depth': 300, 'mutant_tally': 3, 'allele_fraction': 0, 'mutant_depth': 16}}},
                '7:1-5': {'site_depth': 100, 'total_depth': 500, 'Name': 'MUT-BIG>AVE', 'mutant_tally': 4, 'total_mutant_depth': 120, 'total_sites': 5, 
                          'indels': {1: {'site_depth': 300, 'mutant_tally': 3, 'allele_fraction': 0, 'mutant_depth': 110}, -1: {'site_depth': 100, 'mutant_tally': 1, 'allele_fraction': 0, 'mutant_depth': 10}}}, 
               '7:7-11': {'site_depth': 0, 'total_depth': 0, 'Name': 'NO-COV', 'mutant_tally': 0, 'total_mutant_depth': 0, 'total_sites': 5, 'indels': {}}, 
               '8:1-5': {'site_depth': 50, 'total_depth': 250, 'Name': 'MUT-BIG<AVE', 'mutant_tally': 1, 'total_mutant_depth': 49, 'total_sites': 5, 
                          'indels': {-1: {'site_depth': 50, 'mutant_tally': 1, 'allele_fraction': 0, 'mutant_depth': 49}}}, 
               }
#'1:1-5' == wt
#'1:7-11' == wt biggest peak
#'7:1-5' == mut biggest peak, wt_tally != total_sites
#'7:7-11' == no coverage
#'8:1-5' mut buggest peak, but mutant depth < average depth
OUTPUT= {'1:1-5': {'Standard_Deviation': 0, 'Average_Depth': 100, 'Number_of_Peaks': 1, 'Name': 'WT-ONLY', 'IndelLength:AlleleFraction:SupportingCalls': '0:1.0:100'}, 
         '1:7-11': {'Standard_Deviation': '1.210124', 'Average_Depth': 100, 'Number_of_Peaks': 3, 'Name': 'WT-BIGGEST', 'IndelLength:AlleleFraction:SupportingCalls': '-3:0.228571428571:16 -2:0:0 -1:0:0 0:1.0:70 1:0.2:14'}, 
         '7:1-5': {'Standard_Deviation': '0.552771', 'Average_Depth': 100, 'Number_of_Peaks': 2, 'Name': 'MUT-BIG>AVE', 'IndelLength:AlleleFraction:SupportingCalls': '-1:0.0909090909091:10 0:0.0:0 1:1.0:110'}, 
         '7:7-11': {'Standard_Deviation': 0, 'Average_Depth': 0, 'Number_of_Peaks': 0, 'Name': 'NO-COV', 'IndelLength:AlleleFraction:SupportingCalls': '0:0.0:0'},
         '8:1-5': {'Standard_Deviation': '0.140000', 'Average_Depth': 50, 'Number_of_Peaks': 1, 'Name': 'MUT-BIG<AVE', 'IndelLength:AlleleFraction:SupportingCalls': '-1:1.0:49 0:0.0204081632653:1'}}

class TestFormatter(TestBase):
    """
    Test the msi formatter script
    """
    def testCoords(self):
        
        row =['1', '45795895', '45795905', 'WT-ONLY']
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
        self.maxDiff = None
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
        average_depth1=ceil(float(msi_sites1['total_depth'])/msi_sites1['total_sites'])
        wt_ave1=int(msi_sites1['total_depth']-msi_sites1['total_mutant_depth'])/msi_sites1['total_sites']
        wt_frac1=float(wt_ave1)/average_depth1
        highest_frac1 = analyzer.calc_highest_peak(msi_sites1['indels'], wt_frac1,average_depth1)

        msi_sites2=copy.deepcopy(MSI_SITE_DATA['1:7-11'])
        average_depth2=ceil(float(msi_sites2['total_depth'])/msi_sites2['total_sites'])
        wt_ave2=int(msi_sites2['total_depth']-msi_sites2['total_mutant_depth'])/msi_sites2['total_sites']
        wt_frac2=float(wt_ave2)/average_depth2
        highest_frac2 = analyzer.calc_highest_peak(msi_sites2['indels'], wt_frac2, average_depth2)

        self.assertEqual(1.1000000000000001, highest_frac1)
        self.assertEqual(0.94, highest_frac2)        

    def testCalcNumberPeaks(self):
        """Test that the number of peaks and the peak annotation
        is being calculated/parsed correctly. 
        """
        msi_sites1=copy.deepcopy(MSI_SITE_DATA['7:1-5'])

        average_depth1=ceil(float(msi_sites1['total_depth'])/msi_sites1['total_sites'])
        wt_ave1=int(msi_sites1['total_depth']-msi_sites1['total_mutant_depth'])/msi_sites1['total_sites']
        wt_frac1=float(wt_ave1)/average_depth1
        highest_frac1 = analyzer.calc_highest_peak(msi_sites1['indels'], wt_frac1, average_depth1)
        cutoff=0.05
        peaks = []
        wt_sites=analyzer.calc_wildtype(msi_sites1['indels'].keys(), wt_ave1, wt_frac1, highest_frac1)
        num_peaks, sites=analyzer.calc_number_peaks(msi_sites1['indels'], wt_sites, highest_frac1, cutoff)
        output_peaks=3
        output_site_info={-1: '-1:0.0909090909091:10', 0: '0:0.690909090909:76', 1: '1:1.0:110'}
        self.assertEqual(num_peaks, output_peaks)
        self.assertEqual(sites, output_site_info)

    def testCalcWildType1(self):
        """Test the Wildtype calculations"""
        msi_sites=copy.deepcopy(MSI_SITE_DATA['1:7-11'])
        sites = {}
        average_depth=ceil(float(msi_sites['total_depth'])/msi_sites['total_sites'])
        wt_frac =ceil(float(average_depth-msi_sites['total_mutant_depth'])/average_depth)
        wt_ave=int(average_depth-msi_sites['total_mutant_depth'])
        wt_ave=int(wt_ave)
        sites=analyzer.calc_wildtype(msi_sites['indels'].keys(), wt_ave, wt_frac, wt_frac)
        wt_output={0: '0:1.0:70', -1: '-1:0:0', -3: '-3:0:0', -2: '-2:0:0', 1: '1:0:0'}

        self.assertDictEqual(sites, wt_output)

    def testCalcSTDPeaks(self):
        """Test the standard deviation calculations"""
        peaks=['0:0.863414634146:354', '1:0.0402598525993:17', '-1:0.0855382887727:34', '-2:0.0132135895294:5']
        stdev=analyzer.calc_std_peaks(peaks)
        self.assertEqual(stdev, '0.410894')
    
    def testDefineSites(self):
        """ Test that the sites array is created correctly"""
        set1=analyzer.define_sites([3], {})
        set2=analyzer.define_sites([-3], {})
        set3=analyzer.define_sites([-3,2], {})
        set4=analyzer.define_sites([2,4], {})
        
        
        site_output1={0: '0:0:0', 1: '1:0:0', 2: '2:0:0', 3: '3:0:0'}
        site_output2={-3: '-3:0:0', -2: '-2:0:0', -1: '-1:0:0', 0: '0:0:0'}
        site_output3={-3: '-3:0:0', -2: '-2:0:0', -1: '-1:0:0', 0: '0:0:0', 1: '1:0:0',2: '2:0:0'}
        site_output4={0: '0:0:0', 1: '1:0:0', 2: '2:0:0', 3: '3:0:0', 4: '4:0:0'}
        
        self.assertDictEqual(set1, site_output1)
        self.assertDictEqual(set2, site_output2)
        self.assertDictEqual(set3, site_output3)
        self.assertDictEqual(set4, site_output4)
