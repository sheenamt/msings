#!/bin/bash

BAM=$1;

SAVEPATH=$2
PFX=$3
INTERVALS_FILE=$4;
BEDFILE=$5;
MSI_BASELINE=$6;
REF_GENOME=$7;

mkdir -p $SAVEPATH

echo “Starting Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

echo "Making mpileups" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM | awk '{if($4 != 0) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 

echo "Varscan Readcounts start" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
wait

echo "MSI Analyzer start" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt

echo "MSI calls start" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

echo “Completed Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

exit

#NOTE:
#command to make "top-level" output, run after all individual samples are done.
#msi count_msi_samples $MSI_BASELINE $SAVEPATH -o $SAVEPATH/Combined_MSI.txt



