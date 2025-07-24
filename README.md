# Nextflow-Based NeoantigeN Discovery Pipeline (NeonDisco) at CRMY
![NeoNDisco-logo-v2](docs/assets/NeoNDisco-logo-v2.png)

A highly modular Nextflow-based discovery bioinformatics pipeline for prediction of potentially immunogenic recurrent neoantigens from expanded sources of aberrant processes in tumors using RNA-seq data.

## Table of Contents

- [Background](#background)
- [Features](#features)
- [Installation Notes](#installation-notes)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Modes](#pipeline-modes)
- [Input Requirements](#input-requirements)
- [Output Structure](#output-structure)
- [Configuration](#configuration)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Background

Cancer vaccines are an emerging therapeutic option for tumor diseases with a potential to convert *immune-cold* tumours to *immune-hot* tumours in combination with immune checkpoint inhibitor immunotherapy. They are also of particular interest in the quest for the generation of universal cancer vaccines that can be used  “off-the-shelf", in contrast to other highly personalized approaches such as adoptive cell therapies that are likely to be too resource-intensive and expensive to be scaled up outside first-world countries.

[Nextflow](https://www.nextflow.io/) is a free and open sourced dataflow programming platform that enables designing computational workflows to process and analyze bioinformatics datasets in a massively parallel architecture.

This repository documents the design of an integrated cancer neoantigen discovery pipeline with a focus on **alternative sources of neoantigens beyond SNV/indels**. This pipeline takes some inspiration from the Snakemake-based pipeline that focuses on cancer neoantigen discovery previously published in 2023 (Schäfer et al). The main Nextflow script is written in the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) syntax. The file structure of this repository was soft-cloned from [here](https://github.com/FredHutch/workflow-template-nextflow) but the current file content and scripts have been completely customized and rewritten to fit the goal of this pipeline.

## Features

### Core Capabilities
- **Multi-tool fusion detection**: Integrates Arriba, FusionCatcher, and STAR-Fusion for comprehensive fusion gene identification
- **Dual execution modes**: Supports both personalized and cohort-based neoantigen prediction
- **HLA typing integration**: Automated HLA allotype determination using ARCAS-HLA
- **In silico validation**: Transcript validation using FusionInspector and AGFusion
- **Scalable deployment**: Runs locally or on AWS Batch with automatic resource management
- **Flexible input handling**: Supports both local files and S3-based data processing

### Analysis Workflows
1. **Quality Control & Preprocessing**
   - Read trimming with fastp
   - Two-pass STAR alignment with duplicate marking

2. **Fusion Gene Detection**
   - Multi-algorithm fusion calling
   - Consensus filtering and annotation
   - Recurrent fusion identification across cohorts

3. **Neoantigen Prediction**
   - Sample-level HLA-specific predictions
   - Cohort-level HLA frequency analysis
   - pVACfuse-powered neoepitope prediction

4. **HLA Typing**
   - Standalone HLA typing workflow
   - Integration with neoantigen prediction

## Installation Notes

TO BE WRITTEN

## Quick Start

```bash
# Basic local execution with personalized neoantigen prediction
nextflow run neondisco-main.nf \
  -profile local,personalizedNeo \
  -c config/local.config \
  --inputDir /path/to/fastq/files \
  --inputSource local \
  --outputDir ./results

# AWS Batch execution with shared neoantigen prediction
nextflow run neondisco-main.nf \
  -profile awsbatch,sharedNeo \
  -c config/aws.config \
  --manifestPath manifest.tsv \
  --inputSource s3 \
  --outputDir s3://my-bucket/results
```

## Usage

### Command Syntax
```bash
nextflow run neondisco-main.nf -profile <EXECUTOR,MODE[,RESOURCE]> <--OPTION NAME> <ARGUMENT>
```

### Required Arguments

| Parameter | Description | Example |
|-----------|-------------|---------|
| `-c` | Path to configuration file | `-c config/local.config` |
| `-profile` | Execution and analysis mode | `-profile local,personalizedNeo` |
| `--inputDir` | Local directory with BAM/FASTQ files | `--inputDir /data/samples` |
| `--manifestPath` | Tab-delimited manifest file | `--manifestPath samples.tsv` |
| `--inputSource` | Input source type: `local` or `s3` | `--inputSource local` |

### Optional Arguments

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outputDir` | `./outputs` | Output directory (supports S3 URIs) |
| `--trimReads` | `true` | Enable read trimming on FASTQ input |
| `--hlaTypingOnly` | `false` | Run only HLA typing subworkflow |
| `--includeNeoPepPred` | `true` | Include neopeptide prediction |
| `--deleteIntMedFiles` | `false` | Delete intermediate files after use |
| `--deleteStagedFiles` | `true` (if S3) | Delete staged files after processing |
| `--sampleLevelHLANeoPred` | varies | Sample-level HLA neopeptide prediction |
| `--cohortLevelHLANeoPred` | varies | Cohort-level HLA neopeptide prediction |
| `--recurrentFusionsOnly` | varies | Process only recurrent fusions |
| `--recurrenceThreshold` | `0.005` | Recurrence threshold (0.5%) |

## Pipeline Modes

### Execution Profiles

#### Executors
- **`local`**: Run on local machine or HPC cluster
- **`awsbatch`**: Execute on AWS Batch infrastructure

#### Analysis Modes
- **`personalizedNeo`**: Individual sample-based neoantigen prediction
- **`sharedNeo`**: Cohort-wide shared neoantigen analysis

### Profile Examples
```bash
# Local execution, personalized analysis
-profile local,personalizedNeo

# AWS Batch execution, shared analysis
-profile awsbatch,sharedNeo
```

## Input Requirements

### Input File Types
- **FASTQ files**: Paired-end reads (`*R1*.fastq.gz`, `*R2*.fastq.gz`)
- **BAM files**: Aligned reads with corresponding index files

### Manifest File Format
Required columns for `--manifestPath`:
```
sampleName	read1Path	read2Path	sampleType
Sample001	/path/to/sample001_R1.fastq.gz	/path/to/sample001_R2.fastq.gz	Tumor
Sample002	/path/to/sample002_R1.fastq.gz	/path/to/sample002_R2.fastq.gz	Normal
```

### Input Source Compatibility

| Input Source | Executor | Input Method | Supported |
|--------------|----------|--------------|-----------|
| `local` | `local` | `--inputDir` | ✅ |
| `local` | `local` | `--manifestPath` | ✅ |
| `s3` | `awsbatch` | `--manifestPath` | ✅ |
| `s3` | `local` | `--manifestPath` | ❌ |

## Output Structure

```
outputs/
├── preprocessing/
│   ├── trimmed_reads/
│   └── aligned_bams/
├── fusion_calling/
│   ├── arriba/
│   ├── fusioncatcher/
│   ├── starfusion/
│   └── filtered_fusions/
├── hla_typing/
│   ├── sample_hla_types/
│   └── cohort_hla_summary/
├── neoantigen_prediction/
│   ├── sample_level/
│   ├── cohort_level/
│   └── validated_fusions/
└── reports/
    ├── pipeline_info/
    └── execution_reports/
```

## Configuration

### Local Configuration Example
```nextflow
// config/local.config
process {
    executor = 'local'
    cpus = 8
    memory = '32 GB'
}

params {
    max_cpus = 16
    max_memory = '64 GB'
    max_time = '24.h'
}
```

### AWS Configuration Example
```nextflow
// config/aws.config
process {
    executor = 'awsbatch'
    queue = 'my-batch-queue'
    container = 'my-ecr-repo/neondisco:latest'
}

aws {
    region = 'us-east-1'
    batch {
        cliPath = '/usr/local/bin/aws'
    }
}
```

## Advanced Usage

### HLA Typing Only
```bash
nextflow run neondisco-main.nf \
  -profile local,personalizedNeo \
  -c config/local.config \
  --inputDir /path/to/bam/files \
  --inputSource local \
  --hlaTypingOnly true
```

### Recurrent Fusion Analysis
```bash
nextflow run neondisco-main.nf \
  -profile local,sharedNeo \
  -c config/local.config \
  --manifestPath cohort.tsv \
  --inputSource local \
  --recurrentFusionsOnly true \
  --recurrenceThreshold 0.01
```

### S3 Integration
```bash
nextflow run neondisco-main.nf \
  -profile awsbatch,personalizedNeo \
  -c config/aws.config \
  --manifestPath s3://my-bucket/manifest.tsv \
  --inputSource s3 \
  --outputDir s3://my-bucket/results \
  --deleteStagedFiles true
```

## Troubleshooting

### Common Issues

#### Profile Validation Errors
```
ERROR: Only one executor profile allowed
```
**Solution**: Ensure only one executor (`local` or `awsbatch`) is specified:
```bash
-profile local,personalizedNeo  # ✅ Correct
-profile local,awsbatch,personalizedNeo  # ❌ Invalid
```

#### Input Source Compatibility
```
ERROR: AWS Batch executor is not compatible with local input source
```
**Solution**: Match executor with appropriate input source:
```bash
# Local files with local executor
-profile local,personalizedNeo --inputSource local

# S3 files with AWS Batch executor
-profile awsbatch,personalizedNeo --inputSource s3
```

#### Missing Input Files
```
ERROR: No BAM or FASTQ files found
```
**Solution**: Verify file naming conventions:
- FASTQ: `*R1*.fastq.gz`, `*R2*.fastq.gz`
- BAM: `*.bam` with corresponding `*.bai` files

### Resource Requirements

#### Minimum Requirements
- **CPU**: 8 cores
- **Memory**: 32 GB RAM
- **Storage**: 100 GB temporary space per sample

#### Recommended Requirements
- **CPU**: 16+ cores
- **Memory**: 64+ GB RAM
- **Storage**: 500 GB+ temporary space

### Performance Optimization

#### Local Execution
```nextflow
process {
    withName: 'ALIGN_READS_TWOPASS_STARSAM' {
        cpus = 16
        memory = '64 GB'
    }
    withName: 'CALL_FUSIONS_*' {
        cpus = 8
        memory = '32 GB'
    }
}
```

#### AWS Batch Optimization
```nextflow
process {
    withName: 'ALIGN_READS_TWOPASS_STARSAM' {
        queue = 'high-memory-queue'
        cpus = 16
        memory = '64 GB'
    }
}
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:
- Code style and standards
- Testing procedures
- Pull request process
- Issue reporting

## License

TO BE ADDED LATER
