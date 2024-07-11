#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { INDEX_REFERENCE } from './modules/index_reference'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { FASTQC } from './modules/fastqc'
include { BWA_MEM } from './modules/bwa_mem'
include { SAM_TO_BAM } from './modules/sam_to_bam'
include { IVAR_TRIM } from './modules/ivar_trim'
include { IVAR_VARIANTS } from './modules/ivar_variants'
include { IVAR_CONSENSUS } from './modules/ivar_consensus'
include { MULTIQC } from './modules/multiqc'

// Define input channels
read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*_R{1,2}.fastq.gz")
    .ifEmpty { error "Cannot find any read pairs in ${params.input_dir}" }

// Workflow
workflow {
    INDEX_REFERENCE(params.reference)
    TRIMMOMATIC(read_pairs_ch, params.adapter_file)
    FASTQC(TRIMMOMATIC.out.trimmed_reads)
    BWA_MEM(TRIMMOMATIC.out.trimmed_reads, INDEX_REFERENCE.out.index, params.reference)
    SAM_TO_BAM(BWA_MEM.out.sam)
    IVAR_TRIM(SAM_TO_BAM.out.sorted_bam, params.rsv_scheme)
    IVAR_VARIANTS(IVAR_TRIM.out.trimmed_bam, params.reference)
    IVAR_CONSENSUS(IVAR_TRIM.out.trimmed_bam)

    multiqc_files = Channel.empty()
    multiqc_files = multiqc_files.mix(FASTQC.out.fastqc_results.collect())
    multiqc_files = multiqc_files.mix(TRIMMOMATIC.out.trimmed_reads.collect{it[1]}.flatten())
    multiqc_files = multiqc_files.mix(BWA_MEM.out.log.collect())
    multiqc_files = multiqc_files.mix(SAM_TO_BAM.out.sorted_bam.collect{it[1]}.flatten())
    multiqc_files = multiqc_files.mix(IVAR_TRIM.out.trimmed_bam.collect{it[1]}.flatten())
    multiqc_files = multiqc_files.mix(IVAR_VARIANTS.out.variants.collect())
    multiqc_files = multiqc_files.mix(IVAR_CONSENSUS.out.consensus.collect())

    MULTIQC(multiqc_files.collect())
}