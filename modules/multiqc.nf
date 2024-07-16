process MULTIQC_TRIMMED {
    publishDir "${params.output}/quality_control/multiqc_trimmed", mode: 'copy'

    input:
    path('*_fastqc.zip')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}