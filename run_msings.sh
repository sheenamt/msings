#!/bin/bash

INTERVALS_FILE=$1;
BEDFILE=$2;
MSI_BASELINE=$3;
REF_GENOME=$4;
BAM_LIST=${5:}

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0 
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative     
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2

SAVEPATH=$(dirname $BAM)

for BAM in $BAM_LIST; do
    BAMNAME=$(basename $BAM)
    PFX=${BAMNAME%.}
    
    
    mkdir -p $SAVEPATH/$PFX

    echo “Starting Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

    echo "Making mpileups" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM -l $INTERVALS_FILE | awk '{if($4 != 0) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "Varscan Readcounts start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
    wait

    echo "MSI Analyzer start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFXmsi_run_log.txt;
    
    msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt
    
    echo "MSI calls start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    
    msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX --multiplier=$multiplier msi_min_threshold=$msi_min_threshold msi_max_threshold=$msi_max_threshold -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

    echo “Completed Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

fi

echo "Creating summary analysis file for all samples" >> $SAVEPATH/msi_run_log.txt;
msi count_msi_samples $MSI_BASELINE $SAVEPATH -o $SAVEPATH/Combined_MSI.txt



