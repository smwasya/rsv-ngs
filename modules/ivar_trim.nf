process IVAR_TRIM {
    tag "$sample_id"
    
    publishDir "${params.output_dir}/ivar/ivar_trim", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(bam_file), path(bam_index)
    path rsv_scheme

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.sorted.bam"), path("${sample_id}_trimmed.sorted.bam.bai"), emit: trimmed_bam

    script:
    """
    ivar trim -e -b ${rsv_scheme} -p ${sample_id}_trimmed -i ${bam_file}
    samtools sort -@ $task.cpus -o ${sample_id}_trimmed.sorted.bam ${sample_id}_trimmed.bam
    samtools index ${sample_id}_trimmed.sorted.bam
    """
}