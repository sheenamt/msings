"""
Counts #, percentage, and distribution of indels (not base substitutions) within each range of the bed file

Usage:

msi analyzer /path/to/sample/file /path/to/bed/file -o /path/to/outfile
"""

import os
import csv
import sys
import argparse
import re
import pprint
import natsort

from operator import itemgetter
from itertools import chain, groupby
from numpy import std, array
from collections import defaultdict
from msings.utils import multi_split

def build_parser(parser):
    parser.add_argument('msi_file', 
                        type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to the msi output ')
    parser.add_argument('bedfile', 
                        type=argparse.FileType('rU'),
                        default=sys.stdin,
                        help='Path to tab delimited bed file of format chr start end')
    parser.add_argument('-o', '--outfile', 
                        type=argparse.FileType('w'),
                        help='Name of the output file. msi.txt will be tacked to the end')
    parser.add_argument('-p','--peak_fraction_cutoff',
                        default=0.05,
                        type=float,
                        help='Peak fraction cutoff. Default of 0.05')

def calc_msi_dist(variant,msi_site):
    """
    Compute statistics for each site in sample
    """
    # for variant in msi_sites:
    info={}
    #reference/WT info:
    wildtype_reads=variant['base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads'].split(':')[1]
    msi_site['wildtype_depth']+=int(wildtype_reads)
    #Grab all the deletions
    dels = filter(lambda s: "DEL" in str(s), variant['Misc'])
    #And all the insertions
    ins =  filter(lambda s: "INS" in str(s), variant['Misc'])
    sites=dels + ins

    #Want total depth to reflect all sites, not just DEL/INS
    msi_site['total_depth']+=int(variant['q10_depth'])
    msi_site['site_depth']=int(variant['q10_depth'])
    #Keep tally of sites seen
    msi_site['total_sites']+=1

    #Now process DEL/INS specific info
    #Parse the obvious variant format: 'DEL-3-AAA:1:1:23:1:1:0'
    for entry in sites:
        info=multi_split(entry.strip(), ":-")
        if info[0]=='DEL':
            length=int(info[1])*int(-1)
        elif info[0]=='INS':
            length=int(info[1])
        reads=int(info[3])
        #Want a total mutant depth for this loci 
        msi_site['mutant_depth']=msi_site['mutant_depth']+reads
        #Keep tally of mutants seen 
        msi_site['mutant_tally']+=1
        #If we haven't seen this length of mutation, add it
        try:
            allele_fraction=reads/float(msi_site['site_depth'])
        except ZeroDivisionError:
            allele_fraction=0.0
        if not length in msi_site['indels'].keys():
            msi_site['indels'][length]={'allele_fraction':allele_fraction}
            msi_site['indels'][length]['site_depth']=int(variant['q10_depth'])
            msi_site['indels'][length]['mutant_depth']=reads
            msi_site['indels'][length]['mutant_tally']=1
            #Otherwise, update the counts for this mutation length
        else:   
            msi_site['indels'][length]['mutant_tally']+=1
            msi_site['indels'][length]['site_depth']+=int(variant['q10_depth'])
            msi_site['indels'][length]['mutant_depth']=msi_site['indels'][length]['mutant_depth']+reads
            msi_site['indels'][length]['allele_fraction']=(msi_site['indels'][length]['mutant_depth']/float(msi_site['indels'][length]['site_depth']))
    return msi_site

def calc_wildtype(indels, wildtype_ave_depth, wildtype_fraction, highest_peak):
    """
    """
    mx=int(max(indels))+1
    mn=int(min(indels))
    wildtype=":".join([str(0), str(float(wildtype_fraction)/highest_peak), str(wildtype_ave_depth)])
    sites={0:wildtype}

    for key in range(mn,mx):
        if key not in sites.keys():
            sites[key]=str(key)+':0:0'
    return sites

def calc_highest_peak(info, wt_fraction):
    """
    Calculate the highest peak
    """
    highest=wt_fraction
    for loci, details in info.items():
        if details['allele_fraction'] >= highest:
            highest = details['allele_fraction']
    return highest

def calc_number_peaks(info, sites, highest_peak, cutoff):
    """
    Calculate the number of peaks that are above cutoff
    """
    peaks=0
    for loci, details in info.items():
        #take highest allele fraction and divide each alele fraction by that number
        allele_fraction = details['allele_fraction']/highest_peak
        peak=":".join([str(loci), str(allele_fraction), str(details['mutant_depth'])])
        #Overwrite this site if the allele fraction is more than the cutoff
        sites[loci]=peak
        if details['allele_fraction']>= cutoff:
            peaks+=1

    #wildtype is sites[0], if fraction is above cutoff, it counts as a peak
    if sites[0].split(':')[1] >= cutoff:
        peaks+=1
    return peaks, sites

def calc_std_peaks(peaks):
    """Calculate the standard deviation of the alleles
    """
    std_peaks=[]
    for peak in peaks:
        allele=peak.split(':')
        i=0
        for i in range (0,int(allele[2])):
            std_peaks.append(int(allele[0]))
            i=+1

    a=array(std_peaks)
    stddev=format(std(a),'.6f')
    
    return stddev

def calc_summary_stats(output_info, cutoff):
    """Calculate average read depth, number of peaks, standard
    deviation and report each peak for each msi range in the bed file
    """
    sites={}
    #msi_info is all loci for this chromosome
    for name, info in output_info.items():
        #Set total average depth for this site
        if info['total_depth']!= 0 and info['total_sites'] != 0:
            average_depth=int(info['total_depth']/info['total_sites'])
        else:
            average_depth=0
        #Calculate the wildtype information by comparing mutant depth to average depth
        if average_depth != 0 and info['mutant_depth'] < info['wildtype_depth']:
            wildtype_fraction=float(info['wildtype_depth'])/info['total_depth']
            wildtype_ave_depth=int(info['wildtype_depth'])/info['total_sites']
        else:
            wildtype_fraction, wildtype_depth=0,0
            sites[0]='0:0:0'
        if info['indels']:
            highest_peak = calc_highest_peak(info['indels'], wildtype_fraction)
            sites=calc_wildtype(info['indels'].keys(), wildtype_ave_depth, wildtype_fraction, highest_peak)
            num_peaks, peaks=calc_number_peaks(info['indels'], sites, highest_peak, cutoff)
            stdev=calc_std_peaks(peaks.values())
            #Sort the peak list naturally (-3,-2,-1,0,1,2,3)
            peak_list=(" ").join(str(x) for x in natsort.natsorted(peaks.values()))
        else:
            #if there are no indels, its a wiltype call
            wildtype_fraction=1
            wildtype_depth=info['total_depth']
            sites[0]=(":").join(['0', str(float(wildtype_fraction)), str(wildtype_ave_depth)])
            num_peaks=1
            peak_list=sites[0]
            stdev=0
        output_info[name]={'Name':info['Name'],
                           'Average_Depth':average_depth,
                           'Standard_Deviation':stdev,
                           'Number_of_Peaks':num_peaks,
                           'IndelLength:AlleleFraction:Reads':peak_list}
    return output_info
    
def parse_msi_bedfile(row, msi_sites, output_info):
    msi_loci=row['chrom']+':'+row['start']+"-"+row['end']
    #Create chrom dictionary if not already present
    if not row['chrom'] in msi_sites.keys():
        msi_sites[row['chrom']]={}

    for position in xrange(int(row['start']),int(row['end'])+1):
        msi_sites[row['chrom']][position] = msi_loci

    output_info[msi_loci]={'Name':row['name'],
                           'wildtype_depth':0,
                           'total_depth':0,
                           'total_sites':0,
                           'mutant_depth':0,
                           'mutant_tally':0,
                           'indels':{}}
    
    return msi_sites, output_info

def action(args):
    msifile, bedfile, outfile = args.msi_file,args.bedfile, args.outfile
    # prepare a dictionary of chromosome: {range:set(),msi_range:{'start':'1','end':'20'}}
    # can test for membership efficiently
    cutoff=args.peak_fraction_cutoff
    msi_sites={}
    output_info={}
    output = {}
    for row in csv.DictReader(args.bedfile, delimiter='\t', fieldnames=['chrom','start','end','name']):
        msi_sites, output_info = parse_msi_bedfile(row, msi_sites, output_info)
        
    # now we start processing the sample msi info
    sample_msi=csv.DictReader(args.msi_file, delimiter='\t', restkey='Misc')

    #Evaluate the sample-info vs msi-bed info, grouping on chromosome. 
    for row in sample_msi:
        loci_name = msi_sites[row['chrom']][int(row['position'])]
        output_info[loci_name].update(calc_msi_dist(row, output_info[loci_name]))

    output.update(calc_summary_stats(output_info, cutoff))
    fieldnames=['Position','Name','Average_Depth','Number_of_Peaks','Standard_Deviation','IndelLength:AlleleFraction:Reads']

    writer = csv.DictWriter(args.outfile, fieldnames = fieldnames,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    for key, row in natsort.natsorted(output.iteritems()):
       writer.writerow(dict(row, **{'Position': key}))
    

