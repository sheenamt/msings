#!/bin/bash

# Install gatk to $prefix/bin

set -e

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) version [prefix]"
    exit 1
fi

version=$1
#3.2-2
#3.1-1
#2.4-9
#1.6-13
prefix=$2
srcdir=$prefix/src

if [ -f $prefix/bin/gatk${version}/GenomeAnalysisTK.jar ]
then
    if java -jar $prefix/bin/gatk${version}/GenomeAnalysisTK.jar -h 2>&1 | grep -q ${version}; then
	echo "GATK version $version is already installed in $prefix/bin"
	exit 0
    fi
else 
    cd $srcdir
    wget -N ftp://jester.labmed.uw.edu/src/GenomeAnalysisTK-${version}.tar.bz2 
    tar -xf GenomeAnalysisTK-${version}.tar.bz2
    mkdir -p $prefix/bin/gatk-${version} 
    cp GenomeAnalysisTK.jar $prefix/bin/gatk-${version}
fi

