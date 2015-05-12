#!/bin/bash

set -e

if [[ -z $1 ]]; then
   echo "usage: $(basename $0) abs-path-to-fastqs [branch-or-tag]"
   exit 1
fi

fastqs=$1
branch=$2

data_id=$(basename $fastqs)

outdir='/home/genetics/analysis/'$data_id
echo "fastqs: $fastqs"
echo "output directory: $outdir"

if [[ ! -d $fastqs ]]; then
    echo "Error: the directory $fastqs does not exist"
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

# copy the settings-example.conf and insert the data directory
#sed died with using / as delimter (probably because $fastqs has /), @ as delimiter works
sed "s@datadir =@datadir = $fastqs@" settings-example.conf > settings.conf

echo "now execute --> cd $outdir && scons "
