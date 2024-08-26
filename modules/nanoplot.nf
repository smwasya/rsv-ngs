// modules/nanoplot.nf

process NANOPLOT {
    tag "NanoPlot on $sample_id"
    publishDir "${params.output}/nanoplot", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}"

    script:
    """
    mkdir -p ${sample_id}
    NanoPlot --fastq ${fastq} --outdir ${sample_id}
    """
}