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


def calc_msi_dist(site_info,msi_sites):
    """
    Compute statistics for each site in sample
    """
    # for variant in msi_sites:
    info={}
    for variant in site_info:
        #Grab all the deletions
        dels = filter(lambda s: "DEL" in str(s), variant['Misc'])
        #And all the insertions
        ins =  filter(lambda s: "INS" in str(s), variant['Misc'])
        sites=dels + ins
        for key in msi_sites.keys():
            #If this is an MSI loci, proceed
            if int(variant['position']) in msi_sites[key]['range']:
                #Want total depth to reflect all sites, not just DEL/INS
                msi_sites[key]['total_depth']+=int(variant['q10_depth'])
                msi_sites[key]['site_depth']=int(variant['q10_depth'])
                #Keep tally of sites seen
                msi_sites[key]['total_sites']+=1
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
                    msi_sites[key]['mutant_depth']=msi_sites[key]['mutant_depth']+reads
                    #Keep tally of mutants seen 
                    msi_sites[key]['mutant_tally']+=1
                    #If we haven't seen this length of mutation, add it
                    try:
                        allele_fraction=reads/float(msi_sites[key]['site_depth'])
                    except ZeroDivisionError:
                        allele_fraction=0.0
                    if not length in msi_sites[key]['indels'].keys():
                        msi_sites[key]['indels'][length]={'allele_fraction':allele_fraction}
                        msi_sites[key]['indels'][length]['mutant_depth']=reads
                    #Otherwise, update the counts for this mutation length
                    else:    
                        msi_sites[key]['indels'][length]['allele_fraction']=msi_sites[key]['indels'][length]['allele_fraction']+allele_fraction
                        msi_sites[key]['indels'][length]['mutant_depth']=msi_sites[key]['indels'][length]['mutant_depth']+reads
    return msi_sites

def calc_wildtype(info, wildtype_depth, wildtype_fraction):
    """
    """
    mx=int(max(info.keys()))+1
    mn=int(min(info.keys()))
    wildtype=":".join([str(0), str(wildtype_fraction), str(wildtype_depth)])
    sites={0:wildtype}

    for key in range(mn,mx):
        if key not in sites.keys():
            sites[key]=str(key)+':0:0'
            
    return sites
        
def calc_number_peaks(info, sites, cutoff):
    """
    Calculate the number of peaks that are above cutoff
    """
    peaks=0
    for loci, details in info.items():
        peak=":".join([str(loci), str(details['allele_fraction']), str(details['mutant_depth'])])
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

def calc_summary_stats(msi_sites, cutoff):
    """Calculate average read depth, number of peaks, standard
    deviation and report each peak for each msi range in the bed file
    """
    output={}
    msi_sites=msi_sites
    sites={}
    #msi_info is all loci for this chromosome
    for chrom, msi_info in msi_sites.items():
        #info is the site specific information
        for loci, info in msi_info.items():
            #create position for output file
            position=chrom+':'+info['start']+'-'+info['end']
            #Set total average depth for this site
            if info['total_depth']!= 0 and info['total_sites'] != 0:
                average_depth=int(info['total_depth']/info['total_sites'])
            else:
                average_depth=0
            #Calculate the wildtype information by comparing mutant depth to average depth
            if average_depth != 0 and info['mutant_depth'] < average_depth:
                wildtype_fraction=float(average_depth-info['mutant_depth'])/average_depth
                wildtype_depth=int(average_depth-info['mutant_depth'])
            else:
                wildtype_fraction, wildtype_depth=0,0
                sites[0]='0:0:0'
            if info['indels']:
                sites=calc_wildtype(info['indels'], wildtype_depth, wildtype_fraction)
                num_peaks, peaks=calc_number_peaks(info['indels'], sites, cutoff)
                stdev=calc_std_peaks(peaks.values())
                #Sort the peak list naturally (-3,-2,-1,0,1,2,3)
                peak_list=(" ").join(str(x) for x in natsort.natsorted(peaks.values()))
            else:
                #if there are no indels, its a wiltype call
                wildtype_fraction=1
                wildtype_depth=info['total_depth']
                sites[0]=(":").join(['0', str(wildtype_fraction),str(wildtype_depth)])
                num_peaks=1
                peak_list=sites[0]
                stdev=0
            output[position]={'Name':info['name'],
                              'Average_Depth':average_depth,
                              'Total_Depth':info['total_depth'],
                              'Total_Sites':info['total_sites'],
                              'Wildtype_Fraction':wildtype_fraction,
                              'Wildtype_Depth':wildtype_depth,
                              'Standard_Deviation':stdev,
                              'Mutant_Tally':info['mutant_tally'],
                              'Number_of_Peaks':num_peaks,
                              'IndelLength:AlleleFraction:Reads':peak_list}

    return output
    

def parse_msi_bedfile(row, msi_sites):
    msi_loci=row['start']+"-"+row['end']
    #Create chrom dictionary if not already present
    if not row['chrom'] in msi_sites.keys():
        msi_sites[row['chrom']]={}
    #Update msi-loci info for later use
    msi_sites[row['chrom']][msi_loci]={'start':row['start'],
                                       'end':row['end'],
                                       'name':row['name'],
                                       'total_depth':0,
                                       'total_sites':0,
                                       'mutant_depth':0,
                                       'mutant_tally':0,
                                       'indels':{},
                                       'range':set()}
    #Add range of this chromosome for membership testing
    msi_sites[row['chrom']][msi_loci]['range'].update(range(int(row['start']),int(row['end']) + 1))
    return msi_sites
    
def action(args):
    msifile, bedfile, outfile = args.msi_file,args.bedfile, args.outfile
    # prepare a dictionary of chromosome: {range:set(),msi_range:{'start':'1','end':'20'}}
    # can test for membership efficiently
    cutoff=args.peak_fraction_cutoff
    msi_sites={}
    output={}
    for row in csv.DictReader(args.bedfile, delimiter='\t', fieldnames=['chrom','start','end','name']):
        msi_sites.update(parse_msi_bedfile(row, msi_sites))
    # now we start processing the sample msi info
    sample_msi=csv.DictReader(args.msi_file, delimiter='\t', restkey='Misc')
    msi_info = sorted(sample_msi, key=itemgetter('chrom'))
    #Evaluate the sample-info vs msi-bed info, grouping on chromosome. 
    for chrom, site_info in groupby(msi_info, key=itemgetter('chrom')):
        msi_sites[chrom].update(calc_msi_dist(site_info, msi_sites[chrom]))
    output.update(calc_summary_stats(msi_sites, cutoff))

    fieldnames=['Position','Name','Average_Depth','Number_of_Peaks','Standard_Deviation','IndelLength:AlleleFraction:Reads']
    #fieldnames=['Position','Name','Average_Depth','Total_Depth','Total_Sites','Wildtype_Depth','Wildtype_Fraction','Mutant_Tally','Number_of_Peaks','Standard_Deviation','Peak_Distribution']
    #pprint.pprint(output)

    writer = csv.DictWriter(args.outfile, fieldnames = fieldnames,  extrasaction = 'ignore', delimiter = '\t')
    writer.writeheader()
    for key, row in natsort.natsorted(output.iteritems()):
       writer.writerow(dict(row, **{'Position': key}))
    
