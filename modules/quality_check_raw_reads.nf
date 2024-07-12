// Define default parameters
//params.output_dir = 'results'

process QUALITY_CHECK_RAW_READS {
    tag "QC on raw reads: ${sample_id}"
    publishDir "${params.output_dir}/quality_control/fastqc_raw_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*_fastqc.zip"), emit: fastqc_zip
    tuple val(sample_id), path("${sample_id}_*_fastqc.html"), emit: fastqc_html

    script:
    """
    fastqc -o . -f fastq -q ${reads}
    for zip in *_fastqc.zip; do
        mv \$zip ${sample_id}_\$zip
    done
    for html in *_fastqc.html; do
        mv \$html ${sample_id}_\$html
    done
    """
}

process MULTIQC_RAW {
    publishDir "${params.output_dir}/quality_control/multiqc_raw", mode: 'copy'

    input:
    path('*_fastqc.zip')

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data

    script:
    """
    multiqc .
    """
}

workflow QC_RAW_WORKFLOW {
    take:
    reads

    main:
    QUALITY_CHECK_RAW_READS(reads)
    MULTIQC_RAW(QUALITY_CHECK_RAW_READS.out.fastqc_zip.map{ it[1] }.collect())

    emit:
    fastqc_html = QUALITY_CHECK_RAW_READS.out.fastqc_html
    fastqc_zip = QUALITY_CHECK_RAW_READS.out.fastqc_zip
    multiqc_report = MULTIQC_RAW.out.report
    multiqc_data = MULTIQC_RAW.out.data
}
// Example of how to use the workflow in your main script
workflow {
    // Define your input channel here, for example:
    // reads_ch = Channel.fromFilePairs('/path/to/reads/*_{1,2}.fastq.gz')
    
    QC_RAW_WORKFLOW(reads_ch)
}