process EXTRACT_PACBIO_REGIONS {
    
    /*
    This module will extract hg38-mapped pacbio reads based on genomics coordinates specified 
    in the /ref/regions_of_interest_hg38.tsv file
    */
    
    
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