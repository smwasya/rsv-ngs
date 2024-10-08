process ARTIC_MINION {
    tag "$sample_id"
    publishDir "${params.output}/artic_minion/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path(scheme_dir)
    val(scheme)
    val(scheme_version)
    val(medaka_model)

    output:
    tuple val(sample_id), path("${sample_id}.*"), emit: results
    tuple val(sample_id), path("${sample_id}.consensus.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    artic minion \
        --medaka \
        --medaka-model="$params.medaka_model" \
        --normalise 200 \
        --threads 4 \
        --scheme-directory "$params.scheme_dir" \
        --scheme-version 1a RSV_SCHEME \
        --read-file $fastq \
        $sample_id

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic --version 2>&1 | sed 's/artic //')
        medaka: \$(medaka --version 2>&1 | sed 's/medaka //')
    END_VERSIONS
    """
    
}