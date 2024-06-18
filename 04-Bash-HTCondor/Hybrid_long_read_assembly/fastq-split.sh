#!/bin/bash

# arguments from submit file
FASTQ=$1
SPLIT_FASTQ=$2
NUMBER_READS_PER_SPLIT=$3

# determine basename of FASTQ to use as split prefix
FASTQ_BASENAME=`basename $FASTQ .fastq.gz`

# install pigz
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda \
&& miniconda/bin/conda install -y -c conda-forge pigz

# set to use miniconda as base python
export PATH=$PWD/miniconda/bin:$PATH
export HOME=$PWD

# copy FASTQ from staging server
cp $FASTQ .

# calculate number of lines to split
# 4 lines per FASTQ file

let "SPLIT_LINES = $NUMBER_READS_PER_SPLIT * 4"

# create folder to hold split files
mkdir split_fastq

# split FASTQ
time pigz -dc $FASTQ_BASENAME.fastq.gz \
| split -l $SPLIT_LINES \
--filter='pigz > split_fastq/$FILE.fastq.gz' \
- $FASTQ_BASENAME

# copy results back to staging server
mkdir -p $SPLIT_FASTQ/
cp split_fastq/* $SPLIT_FASTQ/

# remove input FASTQ
rm $FASTQ_BASENAME.fastq.gz

# remove miniconda
rm Miniconda3-latest-Linux-x86_64.sh