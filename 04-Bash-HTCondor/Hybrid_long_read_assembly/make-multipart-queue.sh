#!/bin/bash

# create a two column text file listing all the files in a directory (source)
# and a destination folder where output files should be copied on staging server
# used to create a multi-job submission queue file
# example invokation:
# bash 26798-make-multipart-queue.sh \
# /staging/groups/oconnor_group/FLYE/STAGE/cy0333.bonito.split/ \
# /staging/groups/oconnor_group/FLYE/STAGE/cy0333.bonito.mapped \
# map-and-extract.txt \
# *f.fastq.gz

# path to folder containing files to process is first arg to script
INPUT_FOLDER=$1

# path to output folder where CHTC processed files should be copied
OUTPUT_FOLDER=$2

# name of file to save on submit server with two-column information
QUEUE_FILE=$3

# filter specific files from folder to include in queue file
INPUT_FILTER=$4

# get list of files in INPUT_FOLDER and concatenate tab and path to output folder
# write to file that will be used in submit file 
ls -1 "$INPUT_FOLDER"/$INPUT_FILTER | xargs -i echo "{},$OUTPUT_FOLDER" > $QUEUE_FILE