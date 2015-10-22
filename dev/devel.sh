#!/bin/bash

# Usage: devel.sh <virtualenv>

set -e

venv=${1-$(basename $(pwd))-env}

# Install dependencies to it
#dev/bootstrap.sh $venv

# Download test data from jester
test_data=/mnt/disk2/com/test_data/tgc_test_data_141029.tar.gz
tarball=$(basename $test_data)

# If test data is available and the tarball is not already downloaded, get it 
# from the fileshare
if [ -f $test_data ]; then 
    if [ ! -f $tarball ]; then
	scp $test_data .
	tar -xf $tarball
	rm $$tarball
	# Create the data.conf file for the test data:
	bams=${tarball%%.*}
	echo $bams
	echo "[specimen_data]">data.conf
	#Grab only read 1 for parsing the name and absolute path from
	for i in `ls -d -1 $bams/*.bam`; do
	    #get the basename of abs path $i, set pfx to the ID before .1. 
	    pfx=$(basename ${i})
	    #Copy path of read 1, sub for read 2
	    #Print this to the data.conf
	    echo "$pfx = $i" >> data.conf 
	done
fi
echo "now cd $outdir and double check the entries in data.conf. 
There cannot be duplicate sample names. 
Test pipeline setup with 'scons -n', before launching a job with 'scons'"
else
    echo "No test data available, please create a data.conf file listing absolute path to bams"
fi



