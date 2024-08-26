#!/usr/bin/R --vanilla

library(optparse)
library(Biostrings)
library(ggplot2)
library(scales)

calculate_atgc_length <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  seq_str <- as.character(sequences[[1]])
  atgc_length <- sum(strsplit(seq_str, "")[[1]] %in% c("A", "T", "G", "C"))
  return(atgc_length)
}

main <- function(fasta_dir, output_plot, plot_width, plot_height) {
  fasta_files <- list.files(fasta_dir, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE)
  
  cat(sprintf("Found %d FASTA files to process\n", length(fasta_files)))
  
  results <- data.frame(
    Sample = character(),
    ATGC_Length = integer(),
    stringsAsFactors = FALSE
  )
  
  for (fasta_file in fasta_files) {
    cat(sprintf("Processing file: %s\n", fasta_file))
    sample_name <- tools::file_path_sans_ext(basename(fasta_file))
    tryCatch({
      atgc_length <- calculate_atgc_length(fasta_file)
      results <- rbind(results, data.frame(
        Sample = sample_name,
        ATGC_Length = atgc_length
      ))
      cat(sprintf("Processed %s: ATGC length = %d\n", sample_name, atgc_length))
    }, error = function(e) {
      cat(sprintf("Error processing %s: %s\n", fasta_file, conditionMessage(e)))
    })
  }
  
  cat(sprintf("Processed %d samples\n", nrow(results)))
  
  # Create the bar plot with blue bars
  plot <- ggplot(results, aes(x = Sample, y = ATGC_Length)) +
    geom_bar(stat = "identity", fill = "blue") +  # Set bar color to blue
    scale_y_continuous(labels = comma_format()) +
    labs(x = "Samples", y = "ATGC Length", title = "ATGC Genome Lengths by Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, max(results$ATGC_Length) * 1.1))  # Add 10% padding to y-axis

  ggsave(output_plot, plot, width = plot_width, height = plot_height)
  cat(sprintf("Saved plot to %s\n", output_plot))
}

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--output_plot"), default = "atgc_genome_lengths.png", help = "Output plot file name"),
  make_option("--width", type = "double", default = 12, help = "Width of the plot in inches"),
  make_option("--height", type = "double", default = 6, help = "Height of the plot in inches")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser, positional_arguments = 1)

fasta_dir <- args$args[1]
output_plot <- args$options$output_plot
plot_width <- args$options$width
plot_height <- args$options$height

main(fasta_dir, output_plot, plot_width, plot_height)