process MOSDEPTH {
    conda 'bioconda::mosdepth'
    tag "Mosdepth on ${sample_id}"
    publishDir "${params.output}/mosdepth", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam_file)

    output:
    tuple val(sample_id), path("${sample_id}*"), emit: mosdepth_output

    script:
    """
    mosdepth -n --fast-mode --by 200 ${sample_id} ${sorted_bam_file}
    """
}