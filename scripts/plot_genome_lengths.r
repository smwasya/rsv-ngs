#!/usr/bin/env Rscript

library(optparse)
library(Biostrings)
library(ggplot2)
library(reshape2)
library(scales)

calculate_lengths <- function(fasta_file) {
  sequences <- readDNAStringSet(fasta_file)
  total_length <- width(sequences)[1]
  seq_str <- as.character(sequences[[1]])
  n_count <- sum(strsplit(seq_str, "")[[1]] == "N")
  atgc_length <- total_length - n_count
  return(c(total_length, atgc_length, n_count))
}

main <- function(fasta_input, output_plot, output_table, plot_width, plot_height) {
  if (file.info(fasta_input)$isdir) {
    fasta_files <- list.files(fasta_input, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE)
  } else {
    fasta_files <- fasta_input
  }
  
  cat(sprintf("Found %d FASTA files to process\n", length(fasta_files)))
  
  results <- data.frame(
    Sample = character(),
    Total_Length = integer(),
    ATGC_Length = integer(),
    N_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (fasta_file in fasta_files) {
    cat(sprintf("Processing file: %s\n", fasta_file))
    sample_name <- tools::file_path_sans_ext(basename(fasta_file))
    tryCatch({
      lengths <- calculate_lengths(fasta_file)
      results <- rbind(results, data.frame(
        Sample = sample_name,
        Total_Length = lengths[1],
        ATGC_Length = lengths[2],
        N_Count = lengths[3]
      ))
      cat(sprintf("Processed %s: Total length = %d, ATGC length = %d, N count = %d\n",
                  sample_name, lengths[1], lengths[2], lengths[3]))
    }, error = function(e) {
      cat(sprintf("Error processing %s: %s\n", fasta_file, conditionMessage(e)))
    })
  }
  
  cat(sprintf("Processed %d samples\n", nrow(results)))
  
  # Reshape data for plotting
  plot_data <- melt(results, id.vars = "Sample", measure.vars = c("Total_Length", "ATGC_Length"))
  
  # Create the bar plot
  plot <- ggplot(plot_data, aes(x = Sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Total_Length" = "blue", "ATGC_Length" = "red"),
                      labels = c("Total Length", "ATGC Length")) +
    scale_y_continuous(labels = comma_format()) +
    labs(x = "Samples", y = "Genome Length", title = "Genome Lengths by Sample",
         fill = "Length Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    coord_cartesian(ylim = c(0, max(plot_data$value) * 1.1))  # Add 10% padding to y-axis

  ggsave(output_plot, plot, width = plot_width, height = plot_height)
  cat(sprintf("Saved plot to %s\n", output_plot))
  
  # Write the table as TSV
  write.table(results, file = output_table, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved table to %s\n", output_table))
}

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--output_plot"), default = "genome_lengths.png", help = "Output plot file name"),
  make_option(c("-t", "--output_table"), default = "genome_stats.tsv", help = "Output table file name"),
  make_option("--width", type = "double", default = 12, help = "Width of the plot in inches"),
  make_option("--height", type = "double", default = 6, help = "Height of the plot in inches")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser, positional_arguments = 1)

fasta_input <- args$args[1]
output_plot <- args$options$output_plot
output_table <- args$options$output_table
plot_width <- args$options$width
plot_height <- args$options$height

main(fasta_input, output_plot, output_table, plot_width, plot_height)