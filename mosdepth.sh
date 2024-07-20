#!/bin/bash

# Check if both arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_bam_files> <output_directory>"
    exit 1
fi

# Assign input arguments to variables
bam_dir="$1"
output_dir="$2"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Run mosdepth for each sorted BAM file
for sorted_bam_file in "$bam_dir"/*_sorted.bam; do
    if [ -f "$sorted_bam_file" ]; then
        filename=$(basename -- "$sorted_bam_file")
        filename_no_ext="${filename%_sorted.bam}"

        echo "Running mosdepth for $filename_no_ext"
        mosdepth "$output_dir/$filename_no_ext" "$sorted_bam_file"
    fi
done