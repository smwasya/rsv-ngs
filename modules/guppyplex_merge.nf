process GUPPYPLEX_MERGE {
    tag "Guppyplex merge on $barcode_dir"
    publishDir "${params.output}/guppyplex_merge", mode: 'copy'

    input:
    path(barcode_dir)

    output:
    tuple val(barcode_dir.baseName), path("${barcode_dir.baseName}_pass.fastq.gz")

    script:
    """


    artic guppyplex --directory ${barcode_dir} --output ${barcode_dir.baseName}_pass.fastq
    gzip ${barcode_dir.baseName}_pass.fastq
    """
}