"""
Test the subcommand scripts
"""

import os
from os import path
import unittest
import logging
import pprint
import csv
import sys
import json
import subprocess
import filecmp

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

MSI_SITE_DATA={ '1': {'1-5': {'end': '5',
                              'indels': {-1: {'allele_fraction': 0.1, 
                                              'mutant_depth': 10}},
                              'mutant_depth': 10,
                              'mutant_tally': 1,
                              'name': 'NAME1',
                              'range': set([1, 2, 3, 4, 5]),
                              'site_depth': 100,
                              'start': '1',
                              'total_depth': 500,
                              'total_sites': 5},
                      '7-11': {'end': '11',
                               'indels': {-3: {'allele_fraction': 1.4,
                                               'mutant_depth': 140},
                                          1: {'allele_fraction': 1.4, 
                                              'mutant_depth': 140}},
                               'mutant_depth': 280,
                               'mutant_tally': 6,
                               'name': 'NAME2',
                               'range': set([7, 8, 9, 10, 11]),
                               'site_depth': 100,
                               'start': '7',
                               'total_depth': 500,
                               'total_sites': 5}},
                '7': {'1-5': {'end': '5',
                              'indels': {-1: {'allele_fraction': 0.2, 
                                              'mutant_depth': 10},
                                         1: {'allele_fraction': 1.7, 
                                             'mutant_depth': 85}},
                              'mutant_depth': 95,
                              'mutant_tally': 4,
                              'name': 'NAME3',
                              'range': set([1, 2, 3, 4, 5]),
                              'site_depth': 50,
                              'start': '1',
                              'total_depth': 250,
                              'total_sites': 5}}}

OUTPUT={ '1:1-5': {'Average_Depth': 100,
                   'Mutant_Tally': 1,
                   'Name': 'NAME1',
                   'Number_of_Peaks': 2,
                   'IndelLength:AlleleFraction:Reads': '-1:0.1:10 0:0.9:90',
                   'Standard_Deviation': '0.300000',
                   'Total_Depth': 500,
                   'Total_Sites': 5,
                   'Wildtype_Depth': 90,
                   'Wildtype_Fraction': 0.9},
         '1:7-11': {'Average_Depth': 100,
                    'Mutant_Tally': 6,
                    'Name': 'NAME2',
                    'Number_of_Peaks': 3,
                    'IndelLength:AlleleFraction:Reads': '-3:1.4:140 -2:0:0 -1:0:0 0:0:0 1:1.4:140',
                    'Standard_Deviation': '2.000000',
                    'Total_Depth': 500,
                    'Total_Sites': 5,
                    'Wildtype_Depth': 0,
                    'Wildtype_Fraction': 0},
         '7:1-5': {'Average_Depth': 50,
                   'Mutant_Tally': 4,
                   'Name': 'NAME3',
                   'Number_of_Peaks': 3,
                   'IndelLength:AlleleFraction:Reads': '-1:0.2:10 0:0:0 1:1.7:85',
                   'Standard_Deviation': '0.613784',
                   'Total_Depth': 250,
                   'Total_Sites': 5,
                   'Wildtype_Depth': 0,
                   'Wildtype_Fraction': 0}}
MSI_LOCI={'1': {'1-5': {'start': '1',
                        'mutant_tally': 0,
                        'total_depth': 0,
                        'end': '5',
                        'name': 'NAME1',
                        'mutant_depth': 0,
                        'total_sites': 0,
                        'range': set([1, 2, 3, 4, 5]),
                        'indels': {}},
                '7-11': {'start': '7',
                         'mutant_tally': 0,
                         'total_depth': 0,
                         'end': '11',
                         'name': 'NAME2',
                         'mutant_depth': 0,
                         'total_sites': 0,
                         'range': set([8, 9, 10, 11, 7]),
                         'indels': {}}},
          '7': {'1-5': {'start': '1',
                        'mutant_tally': 0,
                        'total_depth': 0,
                        'end': '5',
                        'name': 'NAME3',
                        'mutant_depth': 0,
                        'total_sites': 0,
                        'range': set([1, 2, 3, 4, 5]),
                        'indels': {}}}}

                
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

    def testCalcSummaryStats(self):
        """Test MSI summary calculations
        """
        self.maxDiff=None
        output={}
        cutoff=float(0.05)
        output.update(analyzer.calc_summary_stats(MSI_SITE_DATA, cutoff))
        self.assertDictEqual(OUTPUT, output)               

    def testCalcNumberPeaks(self):
        """Test that the number of peaks and the peak annotation
        is being calculated/parsed correctly. 
        """
        msi_sites=dict(MSI_SITE_DATA['1']['1-5'])
        expected_output=dict(OUTPUT['1:1-5'])
        cutoff=[0.05]
        peaks = []
        sites={0:'0:0:0'}
        peaks, sites=analyzer.calc_number_peaks(msi_sites['indels'], sites, cutoff)
        output_peaks=1
        output_site_info={0: '0:0:0', -1: '-1:0.1:10'}
        self.assertEqual(peaks, output_peaks)
        self.assertEqual(sites, output_site_info)

    def testParseMSIBedfile(self):
        """Test that the MSI Bed file is parsed correctly
        """
        msi_loci={}
        for row in csv.DictReader(open(path.join(msi_testfiles, 'test.msi.bed')), delimiter='\t', fieldnames=['chrom','start','end','name']):
            msi_loci.update(analyzer.parse_msi_bedfile(row, msi_loci))
        self.assertDictEqual(msi_loci, MSI_LOCI)
        
    def testCalcMSIDist(self):
        """Test MSI site distribution calculation"""
        msi_dist={}
        for row in csv.DictReader(open(path.join(msi_testfiles, 'test.msi.bed')), delimiter='\t', fieldnames=['chrom','start','end','name']):
            msi_dist.update(analyzer.parse_msi_bedfile(row, msi_dist))
        sample_msi=csv.DictReader(open(path.join(msi_testfiles, 'test.msi_output')), delimiter='\t', restkey='Misc')
        msi_info = sorted(sample_msi, key=itemgetter('chrom'))
        for chrom, site_info in groupby(msi_info, key=itemgetter('chrom')):
            msi_dist[chrom].update(analyzer.calc_msi_dist(site_info, msi_dist[chrom]))
        self.assertDictEqual(msi_dist, MSI_SITE_DATA)

    def testCalcWildType(self):
        """Test the Wildtype calculations"""
        msi_sites=dict(MSI_SITE_DATA)
        sites = {}
        wt_output={-1: '-1:0:0', 0: '0:0:90', 1: '1:0:0'}
        for chrom, msi_info in msi_sites.items():
            for loci, info in msi_info.items():
                average_depth=info['total_depth']/info['total_sites']
                if average_depth != 0 and info['mutant_depth'] < average_depth:
                    wildtype_fraction=float(average_depth-info['mutant_depth'])/average_depth
                    wildtype_depth=int(average_depth-info['mutant_depth'])
                else:
                    wildtype_fraction, wildtype_tally=0,0
                sites=analyzer.calc_wildtype(info['indels'], wildtype_depth, wildtype_fraction)
        self.assertDictEqual(sites, wt_output)               

    def testCalcSTDPeaks(self):
        """Test the standard deviation calculations"""
        peaks=['0:0.863414634146:354', '1:0.0402598525993:17', '-1:0.0855382887727:34', '-2:0.0132135895294:5']
        stdev=analyzer.calc_std_peaks(peaks)
        self.assertEqual(stdev, '0.410894')

    def testCountMSI(self):
        """Test the output creation"""
        expected_msi_output = os.path.join(msi_testfiles, 'expected_test_msi_output')
        created_msi_output = os.path.join(msi_testfiles, 'created_test_msi_output')
        baseline = os.path.join(msi_testfiles, 'testMSIcontrol')
        cmd=["msi", "count_msi_samples", baseline, msi_testfiles, "-o", created_msi_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_msi_output, created_msi_output))
