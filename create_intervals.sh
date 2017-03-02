#!/bin/bash

source msings-env/bin/activate

BEDFILE=$1;

# get basename of bedfile
BEDBASE=
# get path of bedfile
BEDPATH=
# create pfx of interval file
MSI_INTERVALS=$BEDPATH/$BEDBASE.msi_intervals

echo "Sorting MSI bed file"

sort $BEDFILE  >> $BEDPATH/$BEDBASE.sorted.bed

echo "Making MSI intervals file" 

msi formatter $BEDPATH/$BEDBASE.sorted.bed -o $MSI_INTERVALS

echo "Created $MSI_INTERVALS"




