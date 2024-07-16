process IVAR_CONSENSUS {
    tag "$sample_id"
    
    publishDir "${params.output}/ivar/ivar_consensus", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(trimmed_bam), path(trimmed_bam_index)

    output:
    path "${sample_id}_consensus.fa", emit: consensus

    script:
    """
    samtools mpileup -aa -A -d 0 -Q 0 ${trimmed_bam} | \
    ivar consensus -p ${sample_id}_consensus
    """
}