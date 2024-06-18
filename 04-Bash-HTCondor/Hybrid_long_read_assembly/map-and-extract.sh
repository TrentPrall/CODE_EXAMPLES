#!/bin/bash

# map reads to hg38 reference
# extract immune loci of interest
# copy extracted FASTQ.gz files back to staging server

# arguments from submit file
FASTQ=$1
OUTPUT_FOLDER=$2
REF_GENOME=$3
MINIMAP2_PRESET=$4

# determine basename of FASTQ
FASTQ_BASENAME=`basename $FASTQ .fastq.gz`


# copy reference genome for mapping
cp $REF_GENOME .
REF_GENOME_BASENAME=`basename $REF_GENOME`

# copy FASTQ reads
cp $FASTQ .

# map reads to reference
minimap2 -t 6 -L --eqx -ax $MINIMAP2_PRESET \
$REF_GENOME_BASENAME \
`basename $FASTQ` \
| samtools view -Sbt $REF_GENOME_BASENAME \
| samtools sort - -o aln.bam

# index BAM file
samtools index aln.bam






# extract MHC mapped regions (GABBR1-KIFC1)
samtools view -b aln.bam chr6:29,400,000-33,525,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc.fastq.gz

# extract MHC-1A mapped regions (GABBR1-POLR1H)
samtools view -b aln.bam chr6:29,400,000-30,033,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc-1a.fastq.gz

# extract MHC-class I full mapped regions (GABBR1-MCCD1)
samtools view -b aln.bam chr6:29,400,000-31,499,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc-classI.fastq.gz

# extract MHC-E region (POLR1H-POU5F1)
samtools view -b aln.bam chr6:30,033,000-31,130,500 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc-e.fastq.gz

# extract MHC-1B region (POU5F1-MCCD1)
samtools view -b aln.bam chr6:31,130,500-31,499,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc-1b.fastq.gz

# extract MHC-Class II + III region (MCCD1-KIFC1)
samtools view -b aln.bam chr6:31,499,000-33,525,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.mhc-classII.fastq.gz

# extract KIR mapped regions
samtools view -b aln.bam chr19:54,014,000-55,240,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.kir.fastq.gz

# extract FCGR mapped regions
samtools view -b aln.bam chr1:149,580,000-161,880,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.fcgr.fastq.gz

# extract IGH mapped regions
samtools view -b aln.bam chr14:105386437-107079844 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.igh.fastq.gz

# extract IGK mapped regions
samtools view -b aln.bam chr2:88657361-90435368 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.igk.fastq.gz

# extract IGL mapped regions
samtools view -b aln.bam chr11:21826076-23122913 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.igl.fastq.gz

# extract TRA mapped regions
samtools view -b aln.bam chr14:21421904-22752132 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.tra.fastq.gz

# extract TRB mapped regions
samtools view -b aln.bam chr7:142099011-143013287 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.trb.fastq.gz

# extract Natural Killer Complex mapped regions (M6PR - PRB3)
samtools view -b aln.bam chr12:9,090,000-11,425,000 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.nkc.fastq.gz

# Extract M1 reads
samtools view -b aln.bam cy0325_M1:1-5,670,621 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M1.fastq.gz

# Extract M2 reads
samtools view -b aln.bam cy0161_M2:1-5,468,924 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M2.fastq.gz

# Extract M3 reads
samtools view -b aln.bam cy0333_M3:1-5,227,476 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M3.fastq.gz

# Extract M4 reads
samtools view -b aln.bam cy0695_M4:1-4,988,609 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M4.fastq.gz

# Extract M7 reads
samtools view -b aln.bam cy0390_M7:1-5,497,958 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M7.fastq.gz

# Extract M5 reads
samtools view -b aln.bam cy0355_M5_A:1-1,671,918 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M5A1.fastq.gz

samtools view -b aln.bam cy0355_M5_B:1-820,444 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M5B1.fastq.gz

samtools view -b aln.bam cy0355_M5_II:1-1,265,303 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M5II.fastq.gz

samtools view -b aln.bam cy0424_M5_A:1-1,885,705 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M5A2.fastq.gz

samtools view -b aln.bam cy0424_M5_B:1-3,300,084 \
|samtools fastq - \
|reformat.sh qin=33 int=f in=stdin.fq out=$FASTQ_BASENAME.M5B2.fastq.gz

# copy extracted FASTQ to staging server
mkdir -p $OUTPUT_FOLDER/mhc
mkdir -p $OUTPUT_FOLDER/MCM
mkdir -p $OUTPUT_FOLDER/mhc-1a
mkdir -p $OUTPUT_FOLDER/mhc-classI
mkdir -p $OUTPUT_FOLDER/mhc-e
mkdir -p $OUTPUT_FOLDER/mhc-1b
mkdir -p $OUTPUT_FOLDER/mhc-classII
mkdir -p $OUTPUT_FOLDER/kir
mkdir -p $OUTPUT_FOLDER/fcgr
mkdir -p $OUTPUT_FOLDER/igh
mkdir -p $OUTPUT_FOLDER/igk
mkdir -p $OUTPUT_FOLDER/igl
mkdir -p $OUTPUT_FOLDER/tra
mkdir -p $OUTPUT_FOLDER/trb
mkdir -p $OUTPUT_FOLDER/nkc

cp $FASTQ_BASENAME.mhc.fastq.gz $OUTPUT_FOLDER/mhc/
cp $FASTQ_BASENAME.M*.fastq.gz $OUTPUT_FOLDER/mhc/
cp $FASTQ_BASENAME.M*.fastq.gz $OUTPUT_FOLDER/MCM/
cp $FASTQ_BASENAME.mhc-1a.fastq.gz $OUTPUT_FOLDER/mhc-1a
cp $FASTQ_BASENAME.mhc-classI.fastq.gz $OUTPUT_FOLDER/mhc-classI
cp $FASTQ_BASENAME.mhc-e.fastq.gz $OUTPUT_FOLDER/mhc-e
cp $FASTQ_BASENAME.mhc-1b.fastq.gz $OUTPUT_FOLDER/mhc-1b
cp $FASTQ_BASENAME.mhc-classII.fastq.gz $OUTPUT_FOLDER/mhc-classII
cp $FASTQ_BASENAME.kir.fastq.gz $OUTPUT_FOLDER/kir/
cp $FASTQ_BASENAME.fcgr.fastq.gz $OUTPUT_FOLDER/fcgr/
cp $FASTQ_BASENAME.igh.fastq.gz $OUTPUT_FOLDER/igh/
cp $FASTQ_BASENAME.igk.fastq.gz $OUTPUT_FOLDER/igk/
cp $FASTQ_BASENAME.igl.fastq.gz $OUTPUT_FOLDER/igl/
cp $FASTQ_BASENAME.tra.fastq.gz $OUTPUT_FOLDER/tra/
cp $FASTQ_BASENAME.trb.fastq.gz $OUTPUT_FOLDER/trb/
cp $FASTQ_BASENAME.nkc.fastq.gz $OUTPUT_FOLDER/nkc/

# copy BAM file to staging server since mapping is computationally expensive
# cp aln.bam $OUTPUT_FOLDER/$FASTQ_BASENAME.aln.bam
# cp aln.bam.bai $OUTPUT_FOLDER/$FASTQ_BASENAME.aln.bam.bai

# cleanup
rm $REF_GENOME_BASENAME
rm aln.bam
rm aln.bam.bai
rm $FASTQ_BASENAME.*.fastq.gz
rm `basename $FASTQ`