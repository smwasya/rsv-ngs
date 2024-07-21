process MERGE_GENOMES {
    publishDir "${params.output}/ivar/merged_genome", mode: 'copy'

    input:
    path consensus_files
    path merge_script

    output:
    path "merged_genomes.fa", emit: merged_fasta

    script:
    """
    bash ${merge_script} merged_genomes.fa ${consensus_files}
    """
}
