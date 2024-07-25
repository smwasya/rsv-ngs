
 ![nfmainlog](https://github.com/user-attachments/assets/383a505a-690b-42cc-aa98-89dabf854dd3) 

---
<!-- <img src="https://github.com/user-attachments/assets/383a505a-690b-42cc-aa98-89dabf854dd3" width="500" /> -->
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-5914.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


---
## Introduction

**resvgen** is an advanced bioinformatics pipeline designed specifically for Respiratory Syncytial Virus (RSV) genomics analysis. This tool is currently optimized to process paired-end Illumina sequencing data, guiding researchers through the entire workflow from raw reads to consensus sequence generation.

The pipeline is implemented using [Nextflow](https://www.nextflow.io/), ensuring portability and scalability across various computational environments.

#### Features
1. Specialized for RSV genomics.
2. Supports paired-end Illumina short-read sequencing data.
3. Comprehensive workflow from raw reads to consensus.
4. Utilizes Nextflow for portable and scalable execution.

#### Prerequisites
1. [Nextflow](https://www.nextflow.io/) (version 24.04 or later).
2. [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) (for containerized execution)

#### Installation
```bash
git clone https://github.com/your-repo/resvgen.git
cd resvgen

```

####

## Workflow Overview

<!-- ![resvgen_flow_diagram](https://github.com/user-attachments/assets/f6d035e1-6ec3-497b-b6ec-73f7aeffd388) -->
<img src = "https://github.com/user-attachments/assets/f6d035e1-6ec3-497b-b6ec-73f7aeffd388" height = "200" width = "800" />

By default the pipeline currently performs the following;
   
- Sequence Quality Check ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Trimming Primers ([Trimmomatic](https://docs.tinybio.cloud/docs/trimmomatic-tutorial))
- Post trimming Quality Check ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- 














