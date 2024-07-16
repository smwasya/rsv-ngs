process TRIMMOMATIC {
    tag "$sample_id"
    
    publishDir "${params.output}/trimmomatic", mode: 'copy', enabled: params.organize_output, saveAs: { filename ->
        if (filename.endsWith("_trimmed.fastq.gz")) "trimmed/$filename"
        else if (filename.endsWith("_unpaired.fastq.gz")) "unpaired/$filename"
        else filename
    }

    input:
    tuple val(sample_id), path(reads)
    path adapter_file 

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_R{1,2}_unpaired.fastq.gz", emit: unpaired_reads 

    script:
    """
    trimmomatic PE -threads $task.cpus \
    ${reads[0]} ${reads[1]} \
    ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
    ${sample_id}_R2_trimmed.fastq.gz ${sample_id}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    """
}