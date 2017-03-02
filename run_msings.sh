#!/bin/bash

BAM=$1;
SAVEPATH=$2
PFX=$3
INTERVALS_FILE=doc/mSINGS_TCGA.msi_intervals
BEDFILE=doc/mSINGS_TCGA.bed
MSI_BASELINE=$6;

# CHR: REF_GENOME=/mnt/disk2/com/Genomes/hg19_chr/human_g1k_v37.fasta
REF_GENOME=/mnt/disk2/com/Genomes/gatk-bundle/human_g1k_v37.fasta
VARSCAN=msings-env/bin/VarScan.v2.3.7.jar

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0 
#"msi_min_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI negative   
msi_min_threshold=0.2
#"msi_max_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative   
msi_max_threshold=0.2

source msings-env/bin/activate

mkdir -p $SAVEPATH

echo “Starting Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

echo "Making mpileups" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
samtools mpileup -f $REF_GENOME -d 100000 -A -E  $BAM | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX.mpileup 

echo "Varscan Readcounts start" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;
java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $SAVEPATH/$PFX.msi_output &
wait
 
echo "MSI Analyzer start" >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

msi analyzer $SAVEPATH/$PFX.msi_output $BEDFILE -o $SAVEPATH/$PFX.msi.txt

# echo "MSI calls start" >> $SAVEPATH/msi_run_log.txt;
# date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

# msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

echo “Completed Analysis of $PFX” >> $SAVEPATH/msi_run_log.txt;
date +"%D %H:%M" >> $SAVEPATH/msi_run_log.txt;

#msi create_baseline /path/to/my_output -o /path/to/CUSTOM_MSI_BASELINE
exit

#NOTE:
#command to make "top-level" output, run after all individual samples are done.
#msi count_msi_samples $MSI_BASELINE $SAVEPATH -o $SAVEPATH/Combined_MSI.txt



