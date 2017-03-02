#!/bin/bash

INTERVALS_FILE=$1;
BEDFILE=$2;
REF_GENOME=$3;
#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=${4:}

mkdir -p $SAVEPATH

for BAM in $BAM_LIST; do
    SAVEPATH=$(dirname $BAM)
    BAMNAME=$(basename $BAM)
    PFX=${BAMNAME%.}

    echo “Starting Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

    echo "Making mpileups" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM -l $INTERVALS_FILE | awk '{if($4 != 0) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "Varscan Readcounts start" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
    java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
    wait

    echo "MSI Analyzer start" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

    msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt

echo "Creating baseline of all files" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

msi create_baseline

echo “Completed Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

exit

#NOTE:
#command to make "top-level" output, run after all individual samples are done.
#msi count_msi_samples $MSI_BASELINE $SAVEPATH -o $SAVEPATH/Combined_MSI.txt



