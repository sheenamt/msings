"""
Parsers for the top level summary files
pindel, cnv_exon, cnv_gene, snp, msi, quality, clin_flagged
"""
import os
import csv
import sys
import copy
import natsort
import glob

from itertools import ifilter
from operator import itemgetter

from msings import filters
from msings.utils import munge_pfx
import pandas as pd


"""Each function parses a group of sample files for desired information,
grouping based on the variant_keys list,
some include additional annotation headers,
sample counts, and scores calculated based on counts
"""

def parse_msi(files, control_file, specimens, prefixes, variant_keys, multiplier, threshold):
    """Compare the sample-msi output to the baseline file, report
    Total sites, MSI+ sites and msings score"""
    #Grab the MSI files
    msi_files = ifilter(filters.msi_file_finder,files) 
    msi_files=sorted(msi_files)    

    #Grab the MSI Control info
    control_info=pd.read_csv(control_file, delimiter='\t')
    for i in ['unstable_loci', 'passing_loci', 'msi_status', 'msings_score']:
        control_info = control_info.append({'Position': i}, ignore_index=True)

    variant_keys = 'Position'
    #Now parse each MSI file individually
    for pth in msi_files:
        pfx = munge_pfx(pth.fname)
        try:
            mini_pfx=pfx['mini-pfx'] #if genetics format name
        except TypeError:
            mini_pfx=pfx[0] #otherwise its the filename
        prefixes.append(mini_pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            #Calculate whether this is MSI Stable, MSI Unstable, or Not Covered, based
            # on special MSI calculation
            for line in reader:
                #minimum read depth of 30 to be considered as 'covered'
                if int(line['Average_Depth'])>=30:  
                    pos = control_info.loc[control_info['Position']==line['Position']]
                    value = float(pos['Average']) + (multiplier * float(pos['Standard_Deviation']))
                    if int(line['Number_of_Peaks']) >= value:
                        control_info.loc[(control_info['Position']==line['Position']), mini_pfx] = 'Unstable'
                    else:
                        control_info.loc[(control_info['Position']==line['Position']), mini_pfx] = 'Stable'
                else:
                    control_info.loc[(control_info['Position']==line['Position']), mini_pfx] = 'Not Covered'

    #Now that we're done with the control info, let's make a new dataframe with only the info we want
    specimens = control_info.copy(deep=True)
    specimens = specimens.drop(columns=['Standard_Deviation', 'Average', 'Count'])

    #Parse the user defined thresholds:
    if len(threshold) == 1:
        min_thres=float(threshold)
        max_thres=float(threshold)
    elif len(threshold) == 2:
        min_thres=float(min(threshold))
        max_thres=float(max(threshold))
    else:
        sys.exit("Wrong number of -t thresholds given")
        
    #Create the interpretation based on info parsed by comparing to baseline
    for pfx in prefixes:    
        #Determine total loci in this sample
        total_loci = specimens[specimens[pfx]!='Not Covered'].count()[pfx]
        #Determine unstable loci in this sample
        msi_loci = specimens[specimens[pfx]=='Unstable'].count()[pfx]
        #Add this info to the dataframe as an integer
        specimens.loc[(specimens['Position']=='unstable_loci'), pfx] = "{:.0f}".format(msi_loci)
        specimens.loc[(specimens['Position']=='passing_loci'), pfx] = "{:.0f}".format(total_loci)
        
        #Determine the MSI status, based on threshold given at CLI
        try:
            msings_score=float(msi_loci)/total_loci
            if msings_score >= max_thres:
                status = "POS"
            elif msings_score < min_thres:
                status = "NEG"
            # If the score is between thresholds, its indeterminate
            elif min_thres < msings_score < max_thres:
                status = "IND"
        #If no passing loci, its NEG
        except (ZeroDivisionError, TypeError):
            status = "NEG"
            msings_score = None

        #add the status and score to the dataframe
        specimens.loc[(specimens['Position']=='msi_status'), pfx] = status
        specimens.loc[(specimens['Position']=='msings_score'), pfx] = "{0:.4f}".format(msings_score)

    specimens = specimens.fillna(' ')
    return specimens, prefixes, variant_keys            

def parse_total_mutation_burden(specimens, prefixes, files):
    """Filter for counting as total mutation burden, parses the SNP data file and appends info to the MSI specimens dataframe
    Variant_Type:  All coding or splice; exclude intronic, 5'UTR, 3'UTR, intergenic
    UW_Freq  less than 0.005
    1000g_All= -1 (absent)
    ExAC = -1 (absent)
    Var_reads >8
    """
    variant_to_include = ['frameshift insertion',
                          'frameshift deletion',
                          'frameshift block substitution',
                          'stopgain',
                          'stoploss',
                          'nonframeshift insertion',
                          'nonframeshift deletion',
                          'nonframeshift block substitution',
                          'nonsynonymous SNV',
                          'synonymous SNV',
                          'exonic',
                          'splicing']

    snp_files = ifilter(filters.snp_analysis,files) 
    snp_files=sorted(snp_files) 
    
    #Do not die, or add 'tumor_mutation_burden' info if there are no SNP tabs, which will always happen if run outside UW
    if len(snp_files) == 0:
        print "cannot calculate tumor burden outside UW"
        return specimens

    specimens = specimens.append({'Position': 'tumor_mutation_burden'}, ignore_index=True)

    #Determine total mutation burden, for use internally at UW
    for pth in snp_files:
        pfx = munge_pfx(pth.fname)
        try:
            mini_pfx=pfx['mini-pfx'] #if genetics format name
        except TypeError:
            mini_pfx=pfx[0] #otherwise its the filename
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            count = 0
            total_count = 0
            for line in reader:
                total_count+=1
                for vtype in line['Variant_Type'].split(','):
                    if vtype.strip() in variant_to_include and line['1000g_ALL']=='-1' and float(line['UW_Freq'])<=0.005 and line['EXAC']=='-1' and int(line['Var_Reads'])>=8:
                        count+=1
        #TMB is considered as # of SNPs that passed the filter compared to total # of snps reviewed. 
        tmb = str(count)+"/"+str(total_count)
        specimens.loc[(specimens['Position']=='tumor_mutation_burden'), mini_pfx] = tmb

    return specimens
