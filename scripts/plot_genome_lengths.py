#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse

def calculate_lengths(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length = len(record.seq)
        seq_str = str(record.seq).upper()
        n_count = seq_str.count('N')
        atgc_length = total_length - n_count
        return total_length, atgc_length, n_count

def main(fasta_input, output_plot, output_table, plot_width, plot_height):
    samples = []
    total_lengths = []
    atgc_lengths = []
    n_counts = []

    # Handle multiple input files
    fasta_files = []
    for input_path in fasta_input:
        if os.path.isfile(input_path):
            fasta_files.append(input_path)
        elif os.path.isdir(input_path):
            fasta_files.extend([os.path.join(input_path, f) for f in os.listdir(input_path) 
                                if f.endswith((".fasta", ".fa", ".fna"))])
        else:
            print(f"Warning: Input {input_path} is neither a file nor a directory. Skipping.")

    print(f"Found {len(fasta_files)} FASTA files to process")

    # Process each FASTA file
    for fasta_file in fasta_files:
        print(f"Processing file: {fasta_file}")
        sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
        try:
            total_length, atgc_length, n_count = calculate_lengths(fasta_file)
            
            samples.append(sample_name)
            total_lengths.append(total_length)
            atgc_lengths.append(atgc_length)
            n_counts.append(n_count)
            print(f"Processed {sample_name}: Total length = {total_length}, ATGC length = {atgc_length}, N count = {n_count}")
        except Exception as e:
            print(f"Error processing {fasta_file}: {str(e)}")

    print(f"Processed {len(samples)} samples")

    # Create the bar plot
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))

    bar_width = 0.35
    index = range(len(samples))

    plt.bar(index, total_lengths, bar_width, label='Total Length')
    plt.bar([i + bar_width for i in index], atgc_lengths, bar_width, label='ATGC Length')

    plt.xlabel('Samples')
    plt.ylabel('Genome Length')
    plt.title('Genome Lengths by Sample')
    plt.xticks([i + bar_width/2 for i in index], samples, rotation=45, ha='right')
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()

    print(f"Saved plot to {output_plot}")

    # Write the table as TSV
    with open(output_table, 'w') as tsvfile:
        tsvfile.write("Sample\tTotal Length\tATGC Length\tN Count\n")
        for i in range(len(samples)):
            tsvfile.write(f"{samples[i]}\t{total_lengths[i]}\t{atgc_lengths[i]}\t{n_counts[i]}\n")

    print(f"Saved table to {output_table}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot genome lengths and output N count table from FASTA files.")
    parser.add_argument("fasta_input", nargs='+', help="FASTA file(s) or directory containing FASTA files")
    parser.add_argument("-o", "--output_plot", default="genome_lengths.png", help="Output plot file name (default: genome_lengths.png)")
    parser.add_argument("-t", "--output_table", default="genome_stats.tsv", help="Output table file name (default: genome_stats.tsv)")
    parser.add_argument("--width", type=float, default=12, help="Width of the plot in inches (default: 12)")
    parser.add_argument("--height", type=float, default=6, help="Height of the plot in inches (default: 6)")
    args = parser.parse_args()

    main(args.fasta_input, args.output_plot, args.output_table, args.width, args.height)