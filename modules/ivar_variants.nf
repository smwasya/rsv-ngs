process IVAR_VARIANTS {
    tag "$sample_id"
    
    publishDir "${params.output_dir}/ivar_variants", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(trimmed_bam), path(trimmed_bam_index)
    path reference

    output:
    path "${sample_id}.tsv", emit: variants

    script:
    """
    samtools mpileup -aa -A -d 600000 -B -Q 0 ${trimmed_bam} | \
    ivar variants -p ${sample_id} -q 20 -t 0.03 -r ${reference}
    """
}