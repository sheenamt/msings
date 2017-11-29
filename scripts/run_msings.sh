#!/bin/bash

source msings-env/bin/activate
VARSCAN=msings-env/bin/VarScan.v2.3.7.jar

#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=$1
INTERVALS_FILE=$2;
BEDFILE=$3;
REF_GENOME=$4;
MSI_BASELINE=$5;

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0 
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative     
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2


for BAM in `sed '/^$/d' $BAM_LIST`; do
    SAVEPATH=$(dirname $BAM)
    BAMNAME=$(basename $BAM)
    PFX=${BAMNAME%.*}

    mkdir -p $SAVEPATH/$PFX

    echo “Starting Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

    echo "sorting bam of $PFX" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    samtools sort $BAM $SAVEPATH/$PFX/$PFX.sorted && samtools index $SAVEPATH/$PFX/$PFX.sorted.bam

    echo "Making mpileups" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  $SAVEPATH/$PFX/$PFX.sorted.bam -l $INTERVALS_FILE | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "Varscan Readcounts start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX/$PFX.msi_output &
    wait

    echo "MSI Analyzer start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    
    msi analyzer $SAVEPATH/$PFX/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt
    
    echo "MSI calls start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
     
    msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

    echo “Completed Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

done

echo "Creating summary analysis file for all samples" >> $SAVEPATH/msi_run_log.txt;
msi count_msi_samples $MSI_BASELINE $SAVEPATH -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/Combined_MSI.txt



