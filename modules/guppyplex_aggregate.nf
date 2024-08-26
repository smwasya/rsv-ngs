// modules/guppyplex_aggregate.nf

process ARTIC_GUPPYPLEX_1 {
    tag "Guppyplex aggregate on $barcode_dir"
    publishDir "${params.output}/guppyplex_aggregate", mode: 'copy'

    input:
    path(barcode_dir)

    output:
    tuple val(barcode_dir.baseName), path("${barcode_dir.baseName}_aggregated.fastq.gz")

    script:
    """
    artic guppyplex --directory ${barcode_dir} --output ${barcode_dir.baseName}_aggregated.fastq
    gzip ${barcode_dir.baseName}_aggregated.fastq
    """
}