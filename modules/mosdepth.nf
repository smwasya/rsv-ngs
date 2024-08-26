process MOSDEPTH {
    tag "Mosdepth on ${sample_id}"
    publishDir "${params.output}/mosdepth", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index)
    path regions_file

    output:
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"), emit: regions_bed
    tuple val(sample_id), path("${sample_id}*"), emit: mosdepth_output

    script:
    """
    mosdepth \\
        --threads ${params.threads} \\
        --by ${regions_file} \\
        ${sample_id} \\
        ${bam_file}
    """
}