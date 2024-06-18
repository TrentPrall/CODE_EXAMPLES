process splitFastqProcess {
    input:
    path ont_fastq

    output:
    path "${params.splitDir}/*"

    script:
    """
    echo "Processing file: ${ont_fastq}"
    mkdir -p ${params.splitDir}
    seqkit split2 \
        --by-size ${params.bySize} \
        --extension "${params.extension}" \
        --out-dir ${params.splitDir} \
        --threads ${params.threads} \
        ${ont_fastq}

    echo "Contents of ${params.splitDir} after running seqkit:"
    ls -l ${params.splitDir}
    """
}

workflow {
    Channel.fromPath("${launchDir}/reads/*.fastq.gz").set { fastq_files }

    splitFastqProcess(fastq_files)
}
