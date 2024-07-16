process FASTQC_TRIMMED {
    tag "FASTQC_TRIMMED on ${sample_id}"
    publishDir "${params.output}/quality_control/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html

    script:
    """
    fastqc -q ${read1} ${read2}
    """
}