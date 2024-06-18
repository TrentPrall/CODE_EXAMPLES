#!/bin/bash

# run hifiasm on staged ONT/PACBIO reads


# arguments from submit file
ONT_FASTQ=$1
HIFI_FASTQ=$2
OUTPUT_FOLDER=$3
echo "$ONT_FASTQ"
echo "$HIFI_FASTQ"
echo "$OUTPUT_FOLDER"

ONT=$(basename "$ONT_FASTQ")
echo "$ONT"
PACBIO=$(basename "$HIFI_FASTQ")
echo "$PACBIO"

# set prefix
PREFIX=$(echo "$ONT" | cut -d'.' -f1)
PREFIX="$PREFIX.asm"  # Add ".asm" to the prefix


ls -l

# run hifiasm
echo "hifiasm -o $OUTPUT_FOLDER -t6 --ul $ONT $PACBIO" 
hifiasm -o $PREFIX -t6 --ul $ONT $PACBIO

# make output directory in working folder
echo "mkdir $OUTPUT_FOLDER"
mkdir $OUTPUT_FOLDER
# move output files to output folder
echo "mv $PREFIX* $OUTPUT_FOLDER"
mv $PREFIX* $OUTPUT_FOLDER

# convert hifiasm output to fasta scaffolds

awk '/^S/{print ">"$2"\n"$3}' $OUTPUT_FOLDER/$PREFIX.bp.p_ctg.gfa | fold > $OUTPUT_FOLDER.p_contigs.fasta

# compress output folder
echo "tar -czvf $OUTPUT_FOLDER.tar.gz $OUTPUT_FOLDER/"
tar -czvf $OUTPUT_FOLDER.tar.gz $OUTPUT_FOLDER/


# remove excess fastq files 
echo "rm *.fastq.gz"
rm *.fastq.gz

# output is automatically returned to submit node