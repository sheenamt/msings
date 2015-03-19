#!/bin/bash

# Usage: [PYTHON=/path/to/python] bootstrap.sh [virtualenv-name]
#
# Create a virtualenv and install requirements to it
#
# Specify a python interpreter using 'PYTHON=path/to/python boostrap.sh
#
## System dependencies [Try apt-get] ##
# scons           (The pipeline's make engine)
# g++             (c++ compiler)
# cmake 
# libgsl0-dev     (GNU scientific library)
# libncurses5-dev (For curses.h)
##

set -e

if [[ -z $1 ]]; then
    venv=$(basename $(pwd))-env
else
    venv=$1
fi

if [[ -z $PYTHON ]]; then
    PYTHON=$(which python)
fi



SAMTOOLS_VERSION=0.1.18
VARSCAN_VERSION=v2.3.7

# Create a virtualenv using a specified version of the virtualenv
# source.  This also provides setuptools and pip.  

VENV_URL='http://pypi.python.org/packages/source/v/virtualenv'

# Create virtualenv if necessary
if [ ! -f $venv/bin/activate ]; then
    virtualenv --python /usr/local/bin/python ${venv}
    virtualenv --python /usr/local/bin/python --relocatable ${venv}
else
    echo "found existing virtualenv $venv"
fi

source $venv/bin/activate
mkdir -p $venv/src

# full path; set by activate
venv=$VIRTUAL_ENV

echo $venv

# install samtools 
bash dev/install_samtools.sh ${SAMTOOLS_VERSION} $venv

# install VARSCAN_VERSION=v2.3.7
bash dev/install_varscan.sh ${VARSCAN_VERSION} $venv

python setup.py install

$venv/bin/pip install -r requirements.txt
