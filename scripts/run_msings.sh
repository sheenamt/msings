#!/bin/bash
set -e
source msings-env/bin/activate

#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=$1;
BEDFILE=$2;
MSI_BASELINE=$3;
REF_GENOME=$4;

#Check for required variables:
if [ -z "$BAM_LIST" ]; then echo "BAM_LIST is unset" && exit ; else echo "BAM_LIST is set to '$BAM_LIST'"; fi
if [ -z "$BEDFILE" ]; then echo "BEDFILE is unset" && exit ; else echo "BEDFILE is set to '$BEDFILE'"; fi
if [ -z "$MSI_BASELINE" ]; then echo "MSI_BASELINE is unset" && exit ; else echo "MSI_BASELINE is set to '$MSI_BASELINE'"; fi
if [ -z "$REF_GENOME"]; then echo "REF_GENOME is unset" && exit ; else echo "REF_GENOME is set to '$REF_GENOME'"; fi


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
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  -l $BEDFILE $SAVEPATH/$PFX/$PFX.sorted.bam | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "MSI Analyzer start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    
    msi analyzer $SAVEPATH/$PFX/$PFX.mpileup $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt
    
    echo "MSI calls start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
     
    msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

    echo “Completed Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

done

echo "Creating summary analysis file for all samples" >> $SAVEPATH/msi_run_log.txt;
msi count_msi_samples $MSI_BASELINE $SAVEPATH -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/Combined_MSI.txt
 


