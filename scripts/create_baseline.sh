#!/bin/bash

source msings-env/bin/activate
VARSCAN=msings-env/bin/VarScan.v2.3.7.jar
#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=$1
INTERVALS_FILE=$2;
BEDFILE=$3;
REF_GENOME=$4;

for BAM in `sed '/^$/d' $BAM_LIST`; do
    SAVEPATH=$(dirname $BAM)
    BAMNAME=$(basename $BAM)
    PFX=${BAMNAME%.*}

    mkdir -p $SAVEPATH/$PFX

    echo “Starting Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

    echo "Making mpileups" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM -l $INTERVALS_FILE | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "Varscan Readcounts start" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
    java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
    wait

    echo "MSI Analyzer start" >> $SAVEPATH/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

    msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt

done

echo "Creating baseline of all files" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

msi create_baseline $SAVEPATH -o $SAVEPATH/MSI_BASELINE.txt

echo “Completed Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

exit


