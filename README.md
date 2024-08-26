
<!-- ![nfmainlog](https://github.com/user-attachments/assets/383a505a-690b-42cc-aa98-89dabf854dd3) -->

---
<!-- <img src="https://github.com/user-attachments/assets/383a505a-690b-42cc-aa98-89dabf854dd3" width="500" /> -->
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-5914.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


---
## Introduction

**rsv-ngs** is an advanced bioinformatics pipeline designed specifically for Respiratory Syncytial Virus (RSV) genomics analysis. This tool is currently optimized to process paired-end Illumina and ONT sequencing data, guiding researchers through the entire workflow from raw reads to consensus sequence generation.

The pipeline is implemented using [`Nextflow`](https://www.nextflow.io/), ensuring portability and scalability across various computational environments.

#### Features
1. Specialized for RSV genomics.
2. Supports paired-end Illumina short-read sequencing data.
3. Comprehensive workflow from raw reads to consensus.
4. Utilizes Nextflow for portable and scalable execution.

#### Prerequisites
1. [`Nextflow`](https://www.nextflow.io/) (version 24.04 or later).
2. [`Singularity`](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) (for containerized execution)

#### Installation
```bash
git clone https://github.com/smwasya/rsv-ngs.git
cd rsv-ngs

```
##### Running `resvgen` pipeline

```bash
nexflow run rsv-ngs --input /path/to/input_dir/ --output /path/to/output_dir
```
The current options available for resvgen are `--input`  to specify the user `input` directory and `--output` for `output` directory

####

## Workflow Overview

<!-- ![resvgen_flow_diagram](https://github.com/user-attachments/assets/f6d035e1-6ec3-497b-b6ec-73f7aeffd388) -->
<img src = "https://github.com/user-attachments/assets/f6d035e1-6ec3-497b-b6ec-73f7aeffd388" height = "200" width = "800" />

---

By default the pipeline currently performs the following;
   #### Preprocessing
- Sequence Quality Check ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Trimming Adapters ([`Trimmomatic`](https://docs.tinybio.cloud/docs/trimmomatic-tutorial))
- Post trimming Quality Check ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))

 #### Variant Calling
  - Read alignment ([`BWA-MEM`](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html))
  - Amplicon coverage calculation ([`Mosdepth`](https://github.com/brentp/mosdepth)) ;
    Plots using ([`Custom_script.R`](https://github.com/smwasya/resvgen/blob/main/scripts/plot_mosdepth_regions.r))
  - Sort and Index alignments ([`SAMtools`](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html))
  - Trimming Primers  ([`iVar`](https://github.com/andersen-lab/ivar))
  - Calling variants ([`iVar variants`](https://github.com/andersen-lab/ivar))
  - Variants Visualization ([`Custom_SCript.py`](https://github.com/smwasya/resvgen/tree/main/scripts))
  - Consensus genome calling ([`iVar`](https://github.com/andersen-lab/ivar))
  - Genome length Calculation ([`Custom_script.R`](https://github.com/smwasya/resvgen/blob/main/scripts/plot_genome_lengths.r))

Presents QC reports and visualization plots using ([`MultiQC`](http://multiqc.info/))

##### Clade assignment and lineage classification ([`NextClade`](https://github.com/nextstrain/nextclade))

----

#### Quick Start

1. Follow steps outline here to install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) ; Follow the [tutorial](https://singularity-tutorial.github.io/01-installation/)
3. Download the pipeline and run the command:

```bash
nextflow run rsv-ngs --input /path/to/input_dir --output /path/to/output_dir --platform (illumina|ont) 
``` 















