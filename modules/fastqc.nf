process FASTQC {
    tag "$sample_id"
    
    publishDir "${params.output_dir}/fastqc", mode: 'copy', enabled: params.organize_output

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -t $task.cpus ${trimmed_r1} ${trimmed_r2}
    """
}