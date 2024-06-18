process CONTIGS_TO_FASTA {
    
    /*
    This module will convert the gfa output from hifiasm into a fasta containing contig files
    */
    
    tag "${basename}, ${region}"
    label "map_and_extract"
	publishDir "${params.assembly}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

	input:
    tuple path("hifiasm_files/*"), val(basename), val(region)

	output:
    path "${basename}_${region}.p_contigs.fasta"

	shell:
	'''
    awk '/^S/{print ">"$2;print $3}' hifiasm_files/!{basename}_!{region}.bp.p_ctg.gfa \
    | fold > !{basename}_!{region}.p_contigs.fasta
	'''
}