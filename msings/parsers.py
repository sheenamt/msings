"""
Parsers for the top level summary files
pindel, cnv_exon, cnv_gene, snp, msi, quality, clin_flagged
"""
import os
import csv
import sys
import copy
import natsort

from itertools import count, groupby, chain, ifilter , izip_longest
from operator import itemgetter

from msings import filters

"""Each function parses a group of sample files for desired information,
grouping based on the variant_keys list,
some include additional annotation headers,
sample counts, and scores calculated based on counts
"""


def parse_msi(files, control_file, specimens, prefixes, variant_keys, multiplier, threshold):
    """Compare the sample-msi output to the baseline file, report
    Total sites, MSI+ sites and msings score"""
    #Grab the MSI Control info
    control_info=csv.DictReader(control_file, delimiter='\t')
    control_info=natsort.natsorted(control_info, key=itemgetter('Position'))
    control_info=[d for d in control_info]

    variant_keys = ['Position',]
    for pth in files:
        pfx = pth.fname.replace('.msi.txt','')
        prefixes.append(pfx)
        with open(os.path.join(pth.dir, pth.fname)) as fname:
            reader = csv.DictReader(fname, delimiter='\t')
            sample_msi = natsort.natsorted(reader, key=itemgetter('Position'))
            for key, group in groupby(sample_msi, key=itemgetter('Position')):
                control_row=[d for d in control_info if d['Position']==key]
                try:
                    variant = tuple(control_row[0][k] for k in variant_keys)    
                except IndexError:
                    continue
                for sample_info in group:
                    if int(sample_info['Average_Depth']) >= 30:
                        value = float(control_row[0]['Average_Number_Peaks']) + (multiplier * float(control_row[0]['Standard_Deviation']))
                        if int(sample_info['Number_of_Peaks']) >= value:
                            new_info = 1
                        else:
                            new_info = 0
                    else:           
                        new_info = None
                    specimens[variant][pfx] = new_info
    #Make copy of dictionary to iterate through
    info=copy.copy(specimens)

    #Parse the user defined thresholds:
    if len(threshold) == 1:
        min_thres=float(threshold)
        max_thres=float(threshold)
    elif len(threshold) == 2:
        min_thres=float(min(threshold))
        max_thres=float(max(threshold))
    else:
        sys.exit("Wrong number of -t thresholds given")
    for pfx in prefixes:    
        msi_loci=0
        total_loci=0
        loci=tuple(['passing_loci',])
        msi=tuple(['unstable_loci'],)
        score=tuple(['msing_score'],)
        status=tuple(['msi status'],)
        for entry in info.items():
            if entry[1][pfx] is not None:
                total_loci=total_loci + 1
                msi_loci= msi_loci + entry[1][pfx]
        specimens[loci][pfx]=total_loci
        specimens[msi][pfx]=msi_loci
        try:
            specimens[score][pfx]="{0:.4f}".format(float(msi_loci)/total_loci)
            #POS if above or equal to max threshold
            if float(specimens[score][pfx]) >= max_thres:
                specimens[status][pfx]="POS"
            elif float(specimens[score][pfx]) < min_thres:
                specimens[status][pfx]="NEG"
            # If the score is between thresholds, its indeterminate
            elif min_thres < float(specimens[score][pfx]) < max_thres:
                specimens[status][pfx]="IND"
        except ZeroDivisionError:
            specimens[status][pfx]="NEG"
    fieldnames = variant_keys + list(prefixes) 

    return specimens, prefixes, fieldnames, variant_keys            

