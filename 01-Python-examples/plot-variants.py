#!/opt/anaconda3/envs/genomic-tools/bin/python3

import os, csv, docker

def write_beds(bed):
    ''' takes extracted bed file and creates a single bed file for each entry '''
    tmp = os.getcwd() + '/tmp/'

    ### --- WRITE INVIDUAL BED PER SINGLE ENTRY IN BED FILE --- ###
    with open(bed, 'r') as bed_object:
        bed_reader = csv.reader(bed_object, delimiter='\t')
        for line in bed_reader:
            name = tmp + line[3] + '.tmp.bed'
            with open(name, 'w') as out_object:
                bed_writer = csv.writer(out_object, delimiter='\t')
                bed_writer.writerow(line)


def convert_vcf(bed, vcf):
    ''' takes individual bed files and runs snpsift invtervals and extract fields to prep for karyoploteR '''
    client = docker.from_env() 
    tmp = os.getcwd() + '/tmp/'
    group = vcf.split['.'][-2]
   ### --- loop through entries in bed files --- ###
    with open(bed, 'r') as bed_object:
        bed_reader = csv.reader(bed_object, delimiter='\t')
        for line in bed_reader:
            name = tmp + line[3] + '.tmp.bed'
            vcf_file = tmp + line[3] + '.' + group + '.vcf'
            ### === RUN SNPSIFT INTERVALS --- ###
            client.containers.run('98d1f53855c3', 'intervals ' + bed + ' -i ' + vcf + ' > ' + vcf_file)

            ### --- RUN SNPSIFT EXTRACTFIELDS --- ###
            client.containers.run('98d1f53855c3', 'CHROM POS ID > 23231.baylor12-17.snps.hg38.reorder.casecontrol.goi.clinical.txt'


