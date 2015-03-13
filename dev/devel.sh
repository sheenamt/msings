#!/bin/bash

# Usage: devel.sh <virtualenv>

set -e

venv=${1-$(basename $(pwd))-env}

# Install dependencies to it
dev/bootstrap.sh $venv

# Download test data from jester
wget -N ftp://jester.labmed.uw.edu/genetics/data/tgc_test_data_141029.tar.gz
tar -xf tgc_test_data_141029.tar.gz
rm tgc_test_data_141029.tar.gz

# create settings.conf and use test data
sed 's/datadir =/datadir = 00_TEST_DATA/' settings-example.conf > settings.conf
