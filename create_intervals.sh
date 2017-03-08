#!/bin/bash

source msings-env/bin/activate

BEDFILE=$1;

# get path of bedfile
BEDPATH=$(dirname $BEDFILE)

# get basename of bedfile
BED=$(basename $BEDFILE)
BEDBASE=${BED%.*}

# create pfx of interval file
MSI_INTERVALS=$BEDPATH/$BEDBASE.msi_intervals

echo "Sorting MSI bed file"

sort -V -k1,1 -k2,2n $BEDFILE  >> $BEDPATH/$BEDBASE.sorted.bed

echo "Making MSI intervals file $MSI_INTERVALS" 

msi formatter $BEDPATH/$BEDBASE.sorted.bed $MSI_INTERVALS

echo "Created $MSI_INTERVALS"




