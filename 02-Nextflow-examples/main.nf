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

// Pull modules 
include { EXTRACT_ONT_REGIONS } from "$launchDir/modules/extract_ont_regions.nf"
include { EXTRACT_PACBIO_REGIONS } from "$launchDir/modules/extract_pacbio_regions.nf"
include { MAP_ONT_TO_REF } from "$launchDir/modules/map_ont_to_ref.nf"
include { MAP_PACBIO_TO_REF } from "$launchDir/modules/map_pacbio_to_ref.nf"
include { RUN_HIFIASM } from "$launchDir/modules/run_hifiasm.nf"
include { CONTIGS_TO_FASTA } from "$launchDir/modules/contigs_to_fasta.nf"

// Define Hybrid Assembly workflow
workflow HYBRID_ASSEMBLY {

    // define input channels
    // ont reads supplied by user --ont <absolute path>
    ch_ont_reads = Channel
        .fromPath(params.ont)
        .map { fastq -> tuple(file(fastq), file(fastq).getSimpleName(), "ont") }
        
    // pacbio reads supplied by user --pb <absolute path> 
    ch_pb_reads = Channel
        .fromPath(params.pb)
        .map { fastq -> tuple(file(fastq), file(fastq).getSimpleName(), "pacbio") }
        
    // reference file should be in /ref, path is set in nextflow.config
    ch_ref = Channel.fromPath(params.ref)

    // tsv roi file should be in /ref, path is set in nextflow.config
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
        MAP_ONT_TO_REF(
            ch_ont_reads,
            ch_ref
        )

        // Minimap2 pb hifi to hg38
        MAP_PACBIO_TO_REF(
            ch_pb_reads,
            ch_ref
        )

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

// BEGIN MAIN EXECUTION
workflow {
    HYBRID_ASSEMBLY()
}
