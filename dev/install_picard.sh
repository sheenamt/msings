#!/bin/bash

# Install picard to $prefix/bin

set -e

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) version [prefix]"
    exit 1
fi

version=$1
#1.72
#1.113
prefix=$2
srcdir=$prefix/src

if [ -f $prefix/bin/MarkDuplicates.jar ]
then
    if java -jar $prefix/bin/MarkDuplicates.jar --version 2>&1 | grep -q ${version}; then
	echo "Picard version $version is already installed in $prefix/bin"
	exit 0
    fi
else 
    cd $srcdir
    wget http://sourceforge.net/projects/picard/files/picard-tools/${version}/picard-tools-${version}.zip
    # -o overwrites existing files by the same name and omits interactive prompting
    unzip -j -o picard-tools-${version}.zip -d $prefix/bin

fi
