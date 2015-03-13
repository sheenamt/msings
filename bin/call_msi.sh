#!/bin/bash

#rename variables from commandline
BAM_FILE=$1;
VARSCAN_REF=$2;
VARSCAN=$3;
PFX=$4;
SAVEPATH=$5;
MSI_BED=$6;
MSI_INTERVALS=$7;

echo '==Create mpileup using samtools=='
samtools mpileup -f $VARSCAN_REF -d 1000 -A -E -l $MSI_BED $BAM_FILE | awk '{if ($0 !~ /\t0\t/) {print $0}}' >> $SAVEPATH/$PFX.mpileup

echo '==Clean up in case files already exist=='
rm -f $SAVEPATH/${PFX}.msi_output

echo '==Call counts with varscan=='
java -Xmx4g -jar $VARSCAN readcounts $SAVEPATH/$PFX.mpileup --variants-file $MSI_INTERVALS --min-base-qual 10 --output-file $SAVEPATH/$PFX.msi_output;

echo '==Generate formatted output file=='
perl -I../lib MSI_analyzer.pl $SAVEPATH $PFX $MSI_BED;


