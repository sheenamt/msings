#!/bin/bash

# Install VarScan to $prefix/bin

set -e

if [[ -z $1 ]]; then
    echo "Usage $(basename $0) version [prefix]"
    exit 1
fi

version=$1
#v2.3.6
prefix=$2

if [ -f $prefix/bin/VarScan.${version}.jar ]
then
    echo "Varscan version $version is already installed in $prefix/bin"
    exit 0
else 
    cd $prefix/bin
    wget http://sourceforge.net/projects/varscan/files/VarScan.${version}.jar
    chmod 775 VarScan.${version}.jar
fi

