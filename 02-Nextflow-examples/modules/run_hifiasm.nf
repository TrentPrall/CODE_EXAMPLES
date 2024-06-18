process RUN_HIFIASM {
    
    /*
    This module will run hifiasm using the docker container specified in nextflow.config.
    This is designed to run hybrid assembles on extracted ont and pacbio hifi reads for a single genomic region
    */
    
    
    tag "${basename}, ${region}"
    label "hifiasm"
    publishDir "${params.assembly}/${basename}_${region}", mode: 'copy', overwrite: true

    cpus params.cpus

    input:
        tuple path(pb_fastq), path(ont_fastq), val(basename), val(region)
    
    output:
        tuple path("*"), val(basename), val(region)

        
    script:
	"""
    hifiasm -o ${basename}_${region} -t ${task.cpus} --ul ${ont_fastq} ${pb_fastq}
	"""
}