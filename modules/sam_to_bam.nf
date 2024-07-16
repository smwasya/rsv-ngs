process SAM_TO_BAM {
    tag "$sample_id"
    
    publishDir "${params.output}/sam_to_bam", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: sorted_bam

    script:
    """
    samtools view -@ $task.cpus -bS ${sam_file} | samtools sort -@ $task.cpus -o ${sample_id}_sorted.bam -
    samtools index ${sample_id}_sorted.bam
    """
}