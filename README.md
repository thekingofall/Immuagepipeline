
# ImmuAge Pipeline

Version: 1.0

This pipeline is designed for processing and analyzing RNA-seq data, with a focus on immunological aging studies.

## Prerequisites

- Python 3.x
- R (for certain modules)
- Various bioinformatics tools (STAR, featureCounts, RSEM, etc.)

## Installation

[Provide installation instructions here]

## Usage

```
python immuage_pipeline.py [options]
```

## Options

- `-f, --readFolder`: Folder with RNA-seq data (gzipped fastq files)
- `-e, --threads`: Number of threads for Parafly (Default: 8)
- `-p, --Parallel`: Use Parafly for parallel processing (Default: True)
- `-s, --starindex`: STAR index (Default: '/home/maolp/mao/Ref/hg38/genome')
- `-m, --Module`: Processing module to run (Default: 'RA')
- `-g, --genome`: Path to genome reference (Default: '/home/maolp/mao/Ref/')
- `-u, --gtf`: Path to GTF file (Default: '/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf')
- `-c, --count`: Counting method (Default: 'featureCounts')
- `-a, --alignmethods`: Alignment method (Default: 'hisat2')
- `-t, --alignthreads`: Number of threads for alignment (Default: 40)
- `-x, --rsemindex`: RSEM index (Default: '/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38RSEM/HG38RSEM')
- `-l, --less`: Only process a subset of files (Default: True)
- `-n, --globname`: Global name for the run (Default: 'NEW')
- `-w, --rmarkdown`: Path to R Markdown file (Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuageclock_Run02_PCA.Rmd')
- `-r, --rcut`: Path to R script (Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuage_Rmain.R')
- `--coldata`: Path to coldata file (Required)
- `--add`: Add all, part, or none (Default: 'ap')
- `--health`: Update old program (Default: False)
- `--rootdir`: Root directory (Default: '/data2/users/maolp/Immu_age/Immuageroot')
- `--pipepath`: Pipeline directory (Default: '/home/maolp/mao/Codeman/Immuagepipeline')
- `--resmode`: Result mode (Default: 'train_pre')
- `--config`: Path to config file

## Modules

The pipeline includes several modules that can be specified using the `-m` option:

- S1: FastQC
- S2: TrimQC
- S3_1: STAR indexing
- S3: Alignment
- S14: QC to Count
- S34: Align to count
- S5_1: R Markdown
- S5: Normalization
- S6: Regression
- S78: CircRNA HLA TCR
- S7: TCR HLA
- S8cir: CircleRNA
- RAR: Run all real
- RA: Run all (Default)

## Example

```
python immuage_pipeline.py -f /path/to/fastq/files -m S14 --coldata /path/to/coldata.txt
```

This command will run the pipeline from QC to Count on the fastq files in the specified folder.



