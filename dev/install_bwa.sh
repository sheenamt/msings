#!/bin/bash

# Install bwa to $prefix/bin

set -e

#pipeline v4
#0.6.1

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) version [prefix] "
    exit 1
fi

version=$1

prefix=$2
srcdir=$prefix/src

mkdir -p $srcdir

# Usage: bwa [options]
if $prefix/bin/bwa 2>&1 | grep -q $version
then
    echo "bwa is already installed in $prefix/bin"
    exit 0
else
    cd $srcdir
    wget http://sourceforge.net/projects/bio-bwa/files/bwa-${version}.tar.bz2 
    tar -xjf bwa-$version.tar.bz2 
    make -C bwa-${version}/
    cp bwa-${version}/bwa $prefix/bin
fi   



