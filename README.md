# phasing-pipeline
# SHAPEIT5 Phasing Pipeline for UKBB-like Data

## Table of Contents
1. [Introduction](#introduction)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Input Data Preparation](#input-data-preparation)
5. [Usage](#usage)
6. [Pipeline Steps](#pipeline-steps)
7. [Output](#output)
8. [Troubleshooting](#troubleshooting)
9. [TODOs and Future Improvements](#todos-and-future-improvements)
10. [Contact and Support](#contact-and-support)

## Introduction

This pipeline is designed to perform haplotype phasing on whole-genome sequencing (WGS) data using SHAPEIT5. It is based on the methodology used for the UK Biobank (UKBB) data, including chunking strategies for efficient processing of large datasets and specific handling for the X chromosome.

## Requirements

- Linux operating system (tested on Ubuntu 20.04)
- Bash shell
- sudo privileges for installing dependencies
- At least 32GB RAM (64GB or more recommended for large datasets)
- Sufficient disk space for input, intermediate, and output files (varies with dataset size)

Software dependencies (will be installed by the script):
- SHAPEIT5 (v5.1.1 or later)
- bcftools
- tabix
- wget
- unzip

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your_repo/shapeit5-pipeline.git
   cd shapeit5-pipeline
   ```

2. Make the main script executable:
   ```
   chmod +x shapeit5_pipeline.sh
   ```

3. The script will attempt to install required dependencies. Ensure you have sudo privileges.

## Input Data Preparation

1. Individual VCF files:
   - Place all individual VCF files in a single directory.
   - Files should be gzipped and named consistently (e.g., `sample1.vcf.gz`, `sample2.vcf.gz`, etc.)

2. Genetic map files:
   - Download genetic map files for GRCh38 from the SHAPEIT5 repository.
   - Place these files in a dedicated directory.

3. Pedigree file:
   - Create a pedigree file with three columns: offspring, father, mother.
   - Use 'NA' for unknown parents.

4. Chunking information:
   - Create a file for each chromosome (1-22, X) named `chunks_chr{1..22,X}.txt`.
   - Each file should contain chunk information with start and end positions.
   - Place these files in the same directory as the genetic map files.

5. Haploid males list (for chromosome X):
   - Create a file listing sample IDs of haploid males.
   - Name this file `haploid_males.txt` and place it in an accessible location.

## Usage

Run the pipeline with the following command:

```
./shapeit5_pipeline.sh -i <input_dir> -o <output_dir> -m <map_dir> -p <ped_file> [-t <threads>]
```

Parameters:
- `-i`: Directory containing input VCF files (one per individual)
- `-o`: Directory for output files
- `-m`: Directory containing genetic map files and chunking information
- `-p`: Path to the pedigree file
- `-t`: Number of threads to use (optional, default: 4)

Example:
```
./shapeit5_pipeline.sh -i /path/to/vcfs -o /path/to/output -m /path/to/maps -p /path/to/pedigree.ped -t 8
```

## Pipeline Steps

1. Installation of dependencies
2. Splitting individual VCFs into chromosome-specific files
3. Quality control for each chromosome chunk
4. Phasing of common variants (MAF >= 0.001) in chunks
5. Ligation of phased common variants
6. Phasing of rare variants (MAF < 0.001) in chunks
7. Concatenation of phased rare variants
8. Special processing for chromosome X
9. Cleanup of intermediate files

## Output

The pipeline produces the following main output files:

- `ukb_c{1..22,X}_phased.bcf`: Phased BCF files for each chromosome
- `ukb_c{1..22,X}_phased.bcf.csi`: Index files for the phased BCF files
- `pipeline.log`: Detailed log of the pipeline execution

## Troubleshooting

1. If the script fails to install dependencies, try installing them manually using your system's package manager.

2. For "command not found" errors, ensure that the installed binaries are in your PATH.

3. If you encounter memory errors, try reducing the number of threads or processing chromosomes one at a time.

4. For input/output errors, check disk space and file permissions.

## TODOs and Future Improvements

1. Implement parallelization:
   - Use GNU Parallel or job arrays for cluster environments to process multiple chromosomes simultaneously.

2. Optimize resource management:
   - Implement automatic adjustment of CPU and memory usage based on the system's capabilities.

3. Enhance chromosome X processing:
   - Refine the handling of pseudoautosomal regions (PAR).
   - Improve integration of sex-specific phasing strategies.

4. Add validation steps:
   - Implement quality checks for phased output.
   - Add options for comparing phased data with known pedigrees or reference panels.

5. Improve error handling and reporting:
   - Implement more granular error catching and reporting.
   - Add option for email notifications on pipeline completion or failure.

6. Optimize chunking strategy:
   - Implement adaptive chunking based on variant density and computational resources.

7. Add support for different reference genomes:
   - Extend the pipeline to handle other reference genomes besides GRCh38.

8. Implement checkpointing:
   - Allow the pipeline to resume from the last successful step in case of interruption.

9. Create a configuration file:
   - Move hardcoded parameters to a configurable file for easier customization.

10. Containerization:
    - Create a Docker or Singularity container for easier deployment and reproducibility.

## Contact and Support

For questions, bug reports, or feature requests, please open an issue on the GitHub repository or contact [your_email@example.com].

When reporting issues, please include:
- The exact command you used to run the pipeline
- Relevant sections of the log file
- Any error messages received
- Information about your system (OS, available memory, etc.)
