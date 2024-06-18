#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// log I/O provided by user
log.info        """
                Inputs and Outputs:
                ----------------------------------
                ONT FASTQ              : ${params.ont}
                PacBio FASTQ           : ${params.pb}
                Reference FASTA        : ${params.ref}
                Regions TSV            : ${params.regions_of_interest}
                results_dir            : ${params.results}
                """
                .stripIndent()

workflow HYBRID_ASSEMBLY {

    // define input channels
    ch_ont_reads = Channel
        .fromPath(params.ont)
        .map { fastq -> tuple(file(fastq), file(fastq).getSimpleName(), "ont") }
        

    ch_pb_reads = Channel
        .fromPath(params.pb)
        .map { fastq -> tuple(file(fastq), file(fastq).getSimpleName(), "pacbio") }
        

    ch_ref = Channel.fromPath(params.ref)

    ch_regions_of_interest = Channel
    .fromPath(params.regions_of_interest)
    .splitCsv( header: true, sep: "\t", strip: true )
    .map { 
            row -> tuple( 
                "${row.chromosome}:${row.start}-${row.stop}", row.region, row.merge_key
            ) 
        }

    main:
        // Minimap2 ont to hg38
        MAP_ONT_TO_REF(ch_ont_reads, ch_ref)

        // Minimap2 pb hifi to hg38
        MAP_PACBIO_TO_REF(ch_pb_reads, ch_ref)

        // Extract ONT reads within hg38 genomic intervals 
        EXTRACT_ONT_REGIONS (
            MAP_ONT_TO_REF.out,
            ch_regions_of_interest
        )

        // Extract Pb reads within
        EXTRACT_PACBIO_REGIONS (
            MAP_PACBIO_TO_REF.out,
            ch_regions_of_interest
        )

        // Run hifiasm on extracted reads
        RUN_HIFIASM (
            EXTRACT_ONT_REGIONS.out
            .join ( EXTRACT_PACBIO_REGIONS.out, by: [ 1, 3 ] )
            .map { 
                    basename, region, pb_fastq, pacbio, ont_fastq, ont -> 
                        tuple( file(pb_fastq), file(ont_fastq), basename, region )
                }
        )

        // Convert hifiasm gfa output to fasta of contigs
        CONTIGS_TO_FASTA (
            RUN_HIFIASM.out
        )

}

process MAP_ONT_TO_REF {
    tag "${basename}, ${platform}"
    label "map_and_extract"
    
    cpus params.cpus
    input:
        tuple path(fastq), val(basename), val(platform)
        path reference

    output:
        path "*.bam"

    script:
        minimap2_preset = "map-ont"
        """
        echo "Mapping ${fastq} with ${minimap2_preset} using reference ${reference} for platform ${platform}"
        minimap2 -t ${task.cpus} -L --eqx -ax ${minimap2_preset} \
        ${reference} \
        ${fastq} \
        | samtools view -bS - \
        | samtools sort - -o ${basename}_${platform}.bam
        ls -l ${basename}_${platform}.bam
        """
}

process MAP_PACBIO_TO_REF {
    tag "${basename}, ${platform}"
    label "map_and_extract"
    
    cpus params.cpus
    input:
        tuple path(fastq), val(basename), val(platform)
        path reference

    output:
        path "*.bam"

    script:
        minimap2_preset = "map-hifi"
        """
        echo "Mapping ${fastq} with ${minimap2_preset} using reference ${reference} for platform ${platform}"
        minimap2 -t ${task.cpus} -L --eqx -ax ${minimap2_preset} \
        ${reference} \
        ${fastq} \
        | samtools view -bS - \
        | samtools sort - -o ${basename}_${platform}.bam
        ls -l ${basename}_${platform}.bam
        """
}

process EXTRACT_ONT_REGIONS {
    tag "${basename}, ${platform}, ${region}"
    label "map_and_extract"

    cpus params.cpus

	input:
    each path(bam)
    tuple val(expression), val(region), val(merge_key)

	output:
    tuple path("${basename}_${platform}_${merge_key}.fastq.gz"), val(basename), val(platform), val(merge_key)

	script:
    bam_components = bam.toString().replace(".bam", "").split("_")
    assert bam_components.size() == 2 : "Necessary information could not be parsed from $bam.toString()."
    basename = bam_components[0]
    platform = bam_components[1]
	"""
    samtools index ${bam}
    samtools view -b ${bam} ${expression} \
    | samtools fastq - \
    | reformat.sh qin=33 int=f in=stdin.fq \
    out=${basename}_${platform}_${merge_key}.fastq.gz
	"""
}

process EXTRACT_PACBIO_REGIONS {
    tag "${basename}, ${platform}, ${region}"
    label "map_and_extract"

    cpus params.cpus

	input:
    each path(bam)
    tuple val(expression), val(region), val(merge_key)

	output:
    tuple path("${basename}_${platform}_${merge_key}.fastq.gz"), val(basename), val(platform), val(merge_key)

	script:
    bam_components = bam.toString().replace(".bam", "").split("_")
    assert bam_components.size() == 2 : "Necessary information could not be parsed from $bam.toString()."
    basename = bam_components[0]
    platform = bam_components[1]
	"""
    samtools index ${bam}
    samtools view -b ${bam} ${expression} \
    | samtools fastq - \
    | reformat.sh qin=33 int=f in=stdin.fq \
    out=${basename}_${platform}_${merge_key}.fastq.gz
	"""
}

process RUN_HIFIASM {
    tag "${basename}, ${region}"
    label "hifiasm"
    publishDir "${params.assembly}/${basename}_${region}", mode: 'copy', overwrite: true


    input:
        tuple path(pb_fastq), path(ont_fastq), val(basename), val(region)
    
    output:
        tuple path("*"), val(basename), val(region)

        
    script:
	"""
    hifiasm -o ${basename}_${region} -t ${task.cpus} --ul ${ont_fastq} ${pb_fastq}
	"""
}


process CONTIGS_TO_FASTA {
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
// BEGIN MAIN EXECUTION
workflow {
    HYBRID_ASSEMBLY()
}
