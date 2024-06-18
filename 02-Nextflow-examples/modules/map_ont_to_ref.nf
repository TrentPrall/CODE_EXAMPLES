process MAP_ONT_TO_REF {

    /*
    This module maps ont reads supplied through the -ont command line flag by the user to hg38.
    It sorts the minimap2 sam output and converts the output to bam.
    */


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