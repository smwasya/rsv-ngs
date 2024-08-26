process NEXTCLADE_DATASET {
    

    input:
    val dataset
    val reference
   

    output:
    path "$prefix", emit: dataset

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}"
    def fasta = reference ? "--reference ${reference}" : ''
    
    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        $fasta \\
        $version \\
        --output-dir $prefix
    """
}