process NANOPLOT_TRIMMED {
    tag "NanoPlot on trimmed $sample_id"
    publishDir "${params.output}/nanoplot_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}_trimmed", optional: true

    script:
    """
    if [ \$(zcat ${fastq} | grep -c "^@") -gt 0 ]; then
        mkdir -p ${sample_id}_trimmed
        NanoPlot --fastq ${fastq} --outdir ${sample_id}_trimmed
    else
        echo "No reads found in ${fastq}. Skipping NanoPlot."
        mkdir -p ${sample_id}_trimmed
        echo "No reads found after trimming" > ${sample_id}_trimmed/no_reads.txt
    fi
    """
}