process PLOT_MOSDEPTH_REGIONS {
    label 'process_medium'
    publishDir "${params.output}/mosdepth_plots", mode: 'copy'

    input:
    path(beds)
    path(script)

    output:
    path '*coverage.pdf', emit: coverage_pdf
    path '*coverage.tsv', emit: coverage_tsv
    path '*heatmap.pdf', optional:true, emit: heatmap_pdf
    path '*heatmap.tsv', optional:true, emit: heatmap_tsv
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "mosdepth"
    """
    Rscript $script \
        --input_files ${beds.join(',')} \
        --output_dir . \
        --output_suffix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}