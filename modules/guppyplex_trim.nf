process GUPPYPLEX_TRIM {
    tag "Guppyplex trim on $sample_id"
    publishDir "${params.output}/guppyplex_trim", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_processed.fastq.gz"), emit: processed_fastq
    path "*.log", emit: log

    script:
    """
    echo "Input file:" > ${sample_id}_trim.log
    ls -l ${merged_fastq} >> ${sample_id}_trim.log
    
    echo "Input read count:" >> ${sample_id}_trim.log
    zcat ${merged_fastq} | grep -c "^@" >> ${sample_id}_trim.log
    
    # Create a temporary directory for artic guppyplex
    mkdir -p temp_dir
    cp ${merged_fastq} temp_dir/

    # Run artic guppyplex with error checking
    if ! artic guppyplex --min-length 2500 --max-length 3500 --directory temp_dir --prefix ${sample_id}_trimmed; then
        echo "artic guppyplex failed. Error message:" >> ${sample_id}_trim.log
        artic guppyplex --min-length 2500 --max-length 3500 --directory temp_dir --prefix ${sample_id}_trimmed 2>> ${sample_id}_trim.log
        echo "Using original file due to guppyplex failure." >> ${sample_id}_trim.log
        cp ${merged_fastq} ${sample_id}_processed.fastq.gz
    else
        # Check if any reads passed the filter
        if [ -s "${sample_id}_trimmed.fastq" ]; then
            echo "Reads within 2500-3500 bp found. Using filtered reads." >> ${sample_id}_trim.log
            mv ${sample_id}_trimmed.fastq ${sample_id}_processed.fastq
            gzip ${sample_id}_processed.fastq
        else
            echo "No reads within 2500-3500 bp. Using original file." >> ${sample_id}_trim.log
            cp ${merged_fastq} ${sample_id}_processed.fastq.gz
        fi
    fi

    echo "Output read count:" >> ${sample_id}_trim.log
    zcat ${sample_id}_processed.fastq.gz | grep -c "^@" >> ${sample_id}_trim.log

    # Clean up
    rm -r temp_dir

    # Print directory contents for debugging
    echo "Directory contents:" >> ${sample_id}_trim.log
    ls -l >> ${sample_id}_trim.log
    """
}