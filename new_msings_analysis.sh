#!/bin/bash

set -e

if [[ -z $1 ]]; then
   echo "usage: $(basename $0) abs-path-to-bams [branch-or-tag]"
   exit 1
fi

bams=$1
branch=$2

data_id=$(basename $bams)

outdir='/home/genetics/analysis/'$data_id
echo "bams: $bams"
echo "output directory: $outdir"

if [[ ! -d $bams ]]; then
    echo "Error: the directory $bams does not exist"
    exit 1
fi

#Clone the repo into the data directory
GIT_SSH=/home/genetics/common/ssh/git_ssh.sh \
    git clone git@bitbucket.org:uwlabmed/msings.git $outdir

# git clone does not obey umask setting
chmod g+w $outdir
cd $outdir

# if no branch, tag, or commit is specified, check out the most recent
# tag
if [[ -z $branch ]]; then
   branch=$(git tag | sort -n | tail -1)
fi

echo -n "Using pipeline version $branch"

# Use the specific branch requested
git checkout $branch

# Create the data.conf file for the input data
echo "[specimen_data]">data.conf
#Grab only read 1 for parsing the name and absolute path from
for i in `ls -d -1 $bams/*.bam`; do
    #get the basename of abs path $i, set pfx to the ID before .1. 
    pfx=$(basename ${i})
    #Copy path of read 1, sub for read 2
    #Print this to the data.conf
    echo "$pfx = $i" >> data.conf 
done

echo "now cd $outdir and double check the entries in data.conf. 
There cannot be duplicate sample names. 
Test pipeline setup with 'scons -n', before launching a job with 'scons'"
