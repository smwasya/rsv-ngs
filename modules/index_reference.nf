process INDEX_REFERENCE {
    tag "Indexing reference"
    
    publishDir "${params.output}/index_reference", mode: 'copy', enabled: params.organize_output

    input:
    path reference 

    output:
    path "${reference}.*", emit: index 

    script:
    """
    bwa index ${reference}
    """
}