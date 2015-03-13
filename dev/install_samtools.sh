#!/bin/bash

# Install samtools tools to $prefix/bin

set -e

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) version [prefix]"
    exit 1
fi

version=$1
#0.1.18
prefix=$2
srcdir=$prefix/src

if $prefix/bin/samtools 2>&1 | grep -q ${version}; then
    echo "samtools version $version is already installed in $prefix/bin"
    exit 0
else
    cd $srcdir
    wget -N http://sourceforge.net/projects/samtools/files/samtools/${version}/samtools-${version}.tar.bz2
    tar -xf samtools-${version}.tar.bz2
    make CFLAGS='-g -Wall -O2 -fPIC' -C samtools-${version}
    cp samtools-${version}/samtools $prefix/bin/
fi    



