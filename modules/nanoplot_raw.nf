process NANOPLOT_RAW {
    tag "NanoPlot on raw $sample_id"
    publishDir "${params.output}/nanoplot_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}_raw"

    script:
    """
    mkdir -p ${sample_id}_raw
    NanoPlot --fastq ${fastq} --outdir ${sample_id}_raw
    """
}