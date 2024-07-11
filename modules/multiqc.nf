process MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc .
    """
}