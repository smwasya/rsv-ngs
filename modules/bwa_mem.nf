process BWA_MEM {
    tag "$sample_id"
    
    publishDir "${params.output}/mapping", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)
    path reference_index
    path reference 

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    path "${sample_id}_bwa.log", emit: log 

    script:
    """
    bwa mem -t $task.cpus ${reference} ${trimmed_r1} ${trimmed_r2} > ${sample_id}.sam 2> ${sample_id}_bwa.log
    """
}