#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { INDEX_REFERENCE } from './modules/index_reference'
include { QC_RAW_WORKFLOW } from './modules/quality_check_raw_reads'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { FASTQC_TRIMMED } from './modules/fastqc'
include { BWA_MEM } from './modules/bwa_mem'
include { SAM_TO_BAM } from './modules/sam_to_bam'
include { IVAR_TRIM } from './modules/ivar_trim'
include { IVAR_VARIANTS } from './modules/ivar_variants'
include { IVAR_CONSENSUS } from './modules/ivar_consensus'
include { MULTIQC_TRIMMED } from './modules/multiqc'
include { MOSDEPTH } from './modules/mosdepth'
include { PLOT_MOSDEPTH_REGIONS } from './modules/plot_mosdepth_regions'
include { PLOT_GENOME_LENGTHS } from './modules/plot_genome_lengths'
include { MERGE_GENOMES } from './modules/merge_genomes'

// Validate input parameters
if (!params.input || !params.output) {
    error 'You must specify both --input and --output directories'
}
// Define input channels
Channel
    .fromFilePairs("${params.input}/*_R{1,2}.fastq.gz")
    .ifEmpty { error "Cannot find any read pairs in ${params.input}" }
    .set { read_pairs_ch }

// Define input channels
//read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*_R{1,2}.fastq.gz")
 //   .ifEmpty { error "Cannot find any read pairs in ${params.input_dir}" }

// Workflow
workflow {
    INDEX_REFERENCE(params.reference)
    QC_RAW_WORKFLOW(read_pairs_ch)
    TRIMMOMATIC(read_pairs_ch, params.adapter_file)
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)

    // Run MultiQC on FastQC results
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.zip.map{ it[1] }.collect())

    // Continue with the rest of your workflow
    BWA_MEM(TRIMMOMATIC.out.trimmed_reads, INDEX_REFERENCE.out.index, params.reference)
    SAM_TO_BAM(BWA_MEM.out.sam)
    IVAR_TRIM(SAM_TO_BAM.out.sorted_bam, params.rsv_scheme)
    IVAR_VARIANTS(IVAR_TRIM.out.trimmed_bam, params.reference)
    IVAR_CONSENSUS(IVAR_TRIM.out.trimmed_bam)

    consensus_files = IVAR_CONSENSUS.out.consensus.collect()
    merge_script = file("${projectDir}/scripts/merge_genomes.sh")
    MERGE_GENOMES(consensus_files, merge_script)

    MOSDEPTH(SAM_TO_BAM.out.sorted_bam, file(params.mosdepth_regions))
    mosdepth_regions_files = MOSDEPTH.out.regions_bed.map { it[1] }.collect()
    
    plot_script = file("${projectDir}/scripts/plot_mosdepth_regions.r")
    
    PLOT_MOSDEPTH_REGIONS(mosdepth_regions_files, plot_script)

    
    consensus_dir = IVAR_CONSENSUS.out.consensus.collect().map { it.parent }

    // Run the genome length analysis
    plot_genome_script = file("${projectDir}/scripts/plot_genome_lengths.py")
    PLOT_GENOME_LENGTHS(consensus_dir, plot_genome_script)

    
}