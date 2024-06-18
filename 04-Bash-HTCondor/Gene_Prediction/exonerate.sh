#!/bin/bash

# arguments from submit file
ASSEMBLY=$1
ANIMAL_ID=$2

echo ${ASSEMBLY}
echo ${ANIMAL_ID}

# install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda \
&& miniconda/bin/conda install -c bioconda exonerate \
&& miniconda/bin/conda install -c conda-forge perl
echo 'finished installing exonerate and perl'

export PATH=$PWD/miniconda/bin:$PATH
export HOME=$PWD

# Exonerate with cy0333 annotations
echo 'EXONERATE M3'
exonerate \
--bestn 1 \
--maxintron 20000 \
--showtargetgff \
--showalignment FALSE \
--showvulgar FALSE \
--model affine:local \
--query cy0333_M3_MHC_v2_hybrid_assembly_annotations.fasta \
--target ${ASSEMBLY} \
> ${ANIMAL_ID}_M3_MHC_v2_hybrid_assembly_annotations.gff

echo 'EXONERATE M3 complete'
# process output for exonerate
echo 'EXONERATE M3 process'
bash process-gff.sh \
-e ${ANIMAL_ID}_M3_MHC_v2_hybrid_assembly_annotations.gff \
-p 21295-exonerate_gff_to_alignment_gff3.pl \
-o ${ANIMAL_ID}_M3_MHC_v2_hybrid_assembly_annotations_processed.gff 
echo 'EXONERATE M3 process complete'

# Exonerate with IPD cyno cdna
echo 'EXONERATE cyno IPD'
exonerate \
--bestn 1 \
--maxintron 20000 \
--showtargetgff \
--showalignment FALSE \
--showvulgar FALSE \
--model affine:local \
--query ipd-mhc-mafa-2022-09-19_cleaned.fasta \
--target ${ASSEMBLY} \
> ${ANIMAL_ID}_ipd-mhc-mafa-2022-09-19-cleaned-cdna_annotations.gff
echo 'EXONERATE cyno IPD complete'

# process for geneious import
echo 'EXONERATE cyno IPD process'
bash process-gff.sh \
-e ${ANIMAL_ID}_ipd-mhc-mafa-2022-09-19-cleaned-cdna_annotations.gff \
-p 21295-exonerate_gff_to_alignment_gff3.pl \
-o ${ANIMAL_ID}_ipd-mhc-mafa-2022-09-19-cleaned-cdna_annotations_processed.gff
echo 'EXONERATE cyno IPD process complete'


# Annotate with human HLA CDS
echo 'EXONERATE HLA'
exonerate \
--bestn 1 \
--maxintron 20000 \
--showtargetgff \
--showalignment FALSE \
--showvulgar FALSE \
--model affine:local \
--query hg38_genes_CDS.fasta \
--target ${ASSEMBLY} \
> ${ANIMAL_ID}_hg38_genes_CDS.fasta_annotations.gff
echo 'EXONERATE HLA complete'

# process for geneious import
echo 'EXONERATE HLA process'
bash process-gff.sh \
-e ${ANIMAL_ID}_hg38_genes_CDS.fasta_annotations.gff \
-p 21295-exonerate_gff_to_alignment_gff3.pl \
-o ${ANIMAL_ID}_hg38_genes_CDS.fasta_annotations_processed.gff
echo 'EXONERATE HLA process complete'

# Annotate with Shiina Mamu CDS
echo 'EXONERATE Shiina'
exonerate \
--bestn 1 \
--maxintron 20000 \
--showtargetgff \
--showalignment FALSE \
--showvulgar FALSE \
--model affine:local \
--query Shiina_Mamu_CDS.fasta \
--target ${ASSEMBLY} \
> ${ANIMAL_ID}_Shiina_Mamu_CDS_annotations.gff
echo 'EXONERATE Shiina complete'

# process for geneious import
echo 'EXONERATE Shiina process'
bash process-gff.sh \
-e ${ANIMAL_ID}_Shiina_Mamu_CDS_annotations.gff \
-p 21295-exonerate_gff_to_alignment_gff3.pl \
-o ${ANIMAL_ID}_Shiina_Mamu_CDS_annotations_processed.gff
echo 'EXONERATE Shiina process complete'

# cleanup
echo 'removing excess files'
rm ${ANIMAL_ID}_M3_MHC_v2_hybrid_assembly_annotations.gff
rm ${ANIMAL_ID}_ipd-mhc-mafa-2022-09-19-cleaned-cdna_annotations.gff
rm ${ANIMAL_ID}_hg38_genes_CDS.fasta_annotations.gff
rm ${ANIMAL_ID}_Shiina_Mamu_CDS_annotations.gff

rm *.fasta
rm process-gff.sh
rm 21295-exonerate_gff_to_alignment_gff3.pl
rm Miniconda3-latest-Linux-x86_64.sh

echo 'tar compressing output'
mkdir ${ANIMAL_ID}_MHC_exonerate
mv *.gff ${ANIMAL_ID}_MHC_exonerate
tar -czf ${ANIMAL_ID}_MHC_exonerate.tar.gz ${ANIMAL_ID}_MHC_exonerate/
echo 'job finished sucessfully'
