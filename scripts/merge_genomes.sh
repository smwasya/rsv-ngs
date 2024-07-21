#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_file> <consensus_files...>"
    exit 1
fi

# Extract the output file argument
output_file=$1
shift

# Initialize the output file
: > "$output_file"

# Merge the FASTA files
for fasta in "$@"; do
    cat "$fasta" >> "$output_file"
done
