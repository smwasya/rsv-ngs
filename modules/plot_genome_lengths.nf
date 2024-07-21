process PLOT_GENOME_LENGTHS {
    publishDir "${params.output}/genome_stats", mode: 'copy'

    input:
    path(consensus_dir)
    path(script)

    output:
    path "genome_stats.tsv", emit: stats_tsv
    path "genome_lengths.png", emit: plot_png
    path "debug_output.log", emit: debug_log

    script:
    """
    python3 $script $consensus_dir -o genome_lengths.png -t genome_stats.tsv > debug_output.log 2>&1
    """
}