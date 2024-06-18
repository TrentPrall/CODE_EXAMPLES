#!/bin/bash

# make ONT and PacBio file lists after FASTQ splitting

# ONT
bash make-multipart-queue.sh  \
/staging/groups/oconnor_group/FLYE/STAGE/cy0333/cy0333.bonito.split/ \
/staging/groups/oconnor_group/FLYE/STAGE/cy0333/cy0333.bonito.mapped/ \
ont_list.txt \
*.fastq.gz

# PacBio
bash make-multipart-queue.sh  \
/staging/groups/oconnor_group/FLYE/STAGE/cy0333/cy0333.pacbio.split/ \
/staging/groups/oconnor_group/FLYE/STAGE/cy0333/cy0333.pacbio.mapped/ \
pacbio_list.txt \
*.fastq.gz