#!/bin/bash

# Usage: devel.sh <virtualenv>

set -e

venv=${1-$(basename $(pwd))-env}

# Install dependencies to it
dev/bootstrap.sh $venv

# Download test data from file share
test_data=/mnt/disk2/com/test_data/msings_test_data_141029.tar.gz
tarball=$(basename $test_data)

# If test data is available and the tarball is not already downloaded, get it 
# from the fileshare
if [ -f $test_data ]; then 
    if [ ! -f $tarball ]; then
	scp $test_data .
	tar -xf $tarball
	rm $$tarball
    fi
fi
echo "now cd $outdir and test pipeline setup with scons -n"
