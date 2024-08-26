process COVERAGE_ANALYSIS {
    tag "Coverage analysis"
    publishDir "${params.outdir}/coverage_analysis", mode: 'copy'

    input:
    path(mosdepth_files)

    output:
    path("*.pdf"), emit: pdfs
    path("*.tsv"), emit: tsvs

    script:
    """
    #!/usr/bin/env Rscript

    # Install required packages
    install.packages(c("ggplot2", "scales", "viridis", "tidyverse", "BiocManager"), repos="https://cloud.r-project.org")
    BiocManager::install("ComplexHeatmap")

    # Load required libraries
    library(ggplot2)
    library(scales)
    library(ComplexHeatmap)
    library(viridis)
    library(tidyverse)

    # Function definitions
    escapeRegex <- function(string) {
      gsub("([.|()\\\\^{}+\$*?]|\\\\[|\\\\])", "\\\\\\\\\\\1", string)
    }

    # Set parameters
    INPUT_SUFFIX <- ".regions.bed.gz"
    OUTPUT_SUFFIX <- "regions"
    REGIONS_PREFIX <- NULL

    # Read in data
    INPUT_FILES <- list.files(pattern = paste0(escapeRegex(INPUT_SUFFIX), "\$"), full.names = TRUE)

    if (length(INPUT_FILES) == 0) {
      stop("No *", INPUT_SUFFIX, " files found in the input directory")
    }

    dat <- NULL
    for (input_file in INPUT_FILES) {
      sample <- gsub(INPUT_SUFFIX, '', basename(input_file))
      dat <- rbind(dat, cbind(read.delim(gzfile(input_file), header=FALSE, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)[,-6], sample, stringsAsFactors=F))
    }

    # Reformat table
    if (ncol(dat) == 6) {
      colnames(dat) <- c('chrom', 'start','end', 'region', 'coverage', 'sample')
      if (!is.null(REGIONS_PREFIX)) {
        dat\$region <- as.character(gsub(REGIONS_PREFIX, '', dat\$region))
      }
      dat\$region <- factor(dat\$region, levels=unique(dat\$region[order(dat\$start)]))
    } else {
      colnames(dat) <- c('chrom', 'start','end', 'coverage', 'sample')
    }
    dat\$sample <- factor(dat\$sample, levels=sort(unique(dat\$sample)))

    # Write merged coverage data for all samples to file
    write.table(dat, file="all_samples.regions.coverage.tsv", col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

    # Generate plots
    for (sample in unique(dat\$sample)) {
      sample_dat <- dat[dat\$sample == sample,]
      write.table(sample_dat, file=paste0(sample, ".regions.coverage.tsv"), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
      sample_dat\$coverage <- sample_dat\$coverage + 1

      if (ncol(sample_dat) == 6) {
        plot <- ggplot(sample_dat, aes(x=region, y=coverage)) +
          geom_bar(stat="identity", fill="#D55E00", width=0.6) +
          theme_bw() +
          theme(
            plot.title=element_text(size=10),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=6)) +
          coord_flip() +
          scale_x_discrete(expand=c(0, 0)) +
          scale_y_continuous(
            trans=log10_trans(),
            breaks=10^c(0:10),
            labels=trans_format('log10', math_format(10^.x)),
            expand=c(0, 0)) +
          expand_limits(y=1) +
          ylab(bquote('log'[10]~'(Coverage+1)')) +
          xlab('Amplicon') +
          ggtitle(paste(sample,'median coverage per amplicon'))

        ggsave(file=paste0(sample, ".regions.coverage.pdf"), plot, height=3+(0.2*length(unique(sample_dat\$region))), width=16, units="cm", limitsize=FALSE)
      } else {
        plot <- ggplot(sample_dat, aes(x=end, y=coverage)) +
          geom_ribbon(aes(ymin=0, ymax=coverage), fill="#D55E00") +
          theme_bw() +
          scale_x_continuous(expand=c(0, 0)) +
          scale_y_continuous(
            trans=log10_trans(),
            breaks=10^c(0:10),
            labels=trans_format('log10', math_format(10^.x)),
            expand=c(0, 0)) +
          expand_limits(y=1) +
          ylab(bquote('log'[10]~'(Coverage+1)')) +
          xlab('Position (bp)') +
          ggtitle(paste(sample,'coverage'))

        ggsave(file=paste0(sample, ".regions.coverage.pdf"), plot, height=6, width=12, units="in")
      }
    }

    # Generate heatmap or coverage plot
    if (ncol(dat) == 6) {
      num_samples <- length(unique(dat\$sample))
      
      if (num_samples > 1) {
        mat <- spread(dat[,c("sample", "region", "coverage")], sample, coverage, fill=NA, convert=FALSE)
        rownames(mat) <- mat[,1]
        mat <- as.matrix(mat[,-1])
        mat <- log10(mat + 1)
        
        pdf(file="all_samples.regions.heatmap.pdf", height=10, width=12)
        heatmap(mat, 
                Colv = NA,
                scale = "none",
                main = "Heatmap of log10(Coverage+1) across samples",
                xlab = "Samples",
                ylab = "Regions")
        dev.off()
        
        write.table(cbind(region = rownames(mat), as.data.frame(mat)), 
                    file="all_samples.regions.heatmap.tsv", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
      } else {
        sample_name <- unique(dat\$sample)
        plot <- ggplot(dat, aes(x=region, y=coverage)) +
          geom_bar(stat="identity", fill="#D55E00", width=0.6) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          scale_y_continuous(trans=log10_trans(),
                             breaks=10^c(0:10),
                             labels=trans_format('log10', math_format(10^.x))) +
          labs(title = paste("Coverage across regions for", sample_name),
               x = "Region",
               y = "log10(Coverage+1)")
        
        ggsave(file=paste0(sample_name, ".regions.coverage_plot.pdf"), plot, height=6, width=12)
      }
    }
    """
}