#!/bin/bash

# SHAPEIT5 Phasing Pipeline for UKBB-like Data with Chunking
# This script processes individual VCF files using SHAPEIT5 for phasing, following UKBB methodology

set -e
set -o pipefail

# Function to display script usage
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -m <map_dir> -p <ped_file> [-t <threads>]"
    echo "  -i: Directory containing input VCF files (one per individual)"
    echo "  -o: Directory for output files"
    echo "  -m: Directory containing genetic map files"
    echo "  -p: Pedigree file"
    echo "  -t: Number of threads to use (default: 4)"
    exit 1
}

# Parse command line arguments
while getopts "i:o:m:p:t:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) MAP_DIR="$OPTARG" ;;
        p) PED_FILE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$MAP_DIR" ] || [ -z "$PED_FILE" ]; then
    usage
fi

# Set default number of threads if not provided
THREADS=${THREADS:-4}

# Create output directories
mkdir -p "$OUTPUT_DIR"/{QC,phase_common,ligate,phase_rare,logs,chunks}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$OUTPUT_DIR/logs/pipeline.log"
}

# Function to handle errors
handle_error() {
    log "Error occurred in $1"
    exit 1
}

# Trap errors
trap 'handle_error $LINENO' ERR

# Function to install dependencies
install_dependencies() {
    log "Installing dependencies..."
    sudo apt-get update
    sudo apt-get install -y bcftools tabix wget unzip
    wget https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/shapeit5_v5.1.1_ubuntu20.04_x86_64.zip
    unzip shapeit5_v5.1.1_ubuntu20.04_x86_64.zip
    sudo mv shapeit5_* /usr/local/bin/
    rm shapeit5_v5.1.1_ubuntu20.04_x86_64.zip
    log "Dependencies installed successfully"
}

# Function to split individual VCFs into chromosomes
split_vcfs() {
    log "Splitting individual VCFs into chromosomes"
    for vcf in "$INPUT_DIR"/*.vcf.gz; do
        for chr in {1..22} X; do
            bcftools view -r chr$chr $vcf -Oz -o "$OUTPUT_DIR/chunks/$(basename $vcf .vcf.gz)_chr$chr.vcf.gz"
            bcftools index "$OUTPUT_DIR/chunks/$(basename $vcf .vcf.gz)_chr$chr.vcf.gz"
        done
    done
}

# Function to perform quality control
perform_qc() {
    local chr=$1
    local chunk=$2
    log "Performing QC for chromosome $chr chunk $chunk"
    
    bcftools norm --threads $THREADS -m -any "$OUTPUT_DIR/chunks/ukb_c${chr}_b${chunk}_v1.vcf.gz" -Ou | \
    bcftools annotate -x "FORMAT" -Ou | \
    bcftools +fill-tags --threads $THREADS -Ou -- -t HWE,AF,ExcHet -S "$PED_FILE" | \
    bcftools view --threads $THREADS -f PASS -e 'ALT="*" | F_MISSING>0.1 | INFO/AAScore<0.8 | INFO/AF=0 | INFO/AF=1 | INFO/HWE<1e-30' -Ob -o "$OUTPUT_DIR/QC/ukb_c${chr}_b${chunk}_qc.bcf"
    
    bcftools index "$OUTPUT_DIR/QC/ukb_c${chr}_b${chunk}_qc.bcf"
}

# Function to phase common variants
phase_common() {
    local chr=$1
    local chunk=$2
    log "Phasing common variants for chromosome $chr chunk $chunk"
    
    shapeit5_phase_common --input "$OUTPUT_DIR/QC/ukb_c${chr}_b${chunk}_qc.bcf" \
                          --map "$MAP_DIR/chr${chr}.b38.gmap.gz" \
                          --output "$OUTPUT_DIR/phase_common/ukb_c${chr}_b${chunk}_common.bcf" \
                          --thread $THREADS \
                          --filter-maf 0.001 \
                          --region $chunk \
                          --pedigree "$PED_FILE"
}

# Function to ligate phased common variants
ligate_common() {
    local chr=$1
    log "Ligating phased common variants for chromosome $chr"
    
    ls -1v "$OUTPUT_DIR/phase_common/ukb_c${chr}_b*_common.bcf" > "$OUTPUT_DIR/ligate/list_ligate_chr${chr}.txt"
    
    shapeit5_ligate --input "$OUTPUT_DIR/ligate/list_ligate_chr${chr}.txt" \
                    --output "$OUTPUT_DIR/ligate/ukb_c${chr}_common_ligated.bcf" \
                    --thread $THREADS
}

# Function to phase rare variants
phase_rare() {
    local chr=$1
    local chunk=$2
    log "Phasing rare variants for chromosome $chr chunk $chunk"
    
    shapeit5_phase_rare --input "$OUTPUT_DIR/QC/ukb_c${chr}_b${chunk}_qc.bcf" \
                        --input-scaffold "$OUTPUT_DIR/ligate/ukb_c${chr}_common_ligated.bcf" \
                        --output "$OUTPUT_DIR/phase_rare/ukb_c${chr}_b${chunk}_rare.bcf" \
                        --thread $THREADS \
                        --region $chunk
}

# Function to concatenate phased rare variants
concatenate_rare() {
    local chr=$1
    log "Concatenating phased rare variants for chromosome $chr"
    
    ls -1v "$OUTPUT_DIR/phase_rare/ukb_c${chr}_b*_rare.bcf" > "$OUTPUT_DIR/phase_rare/list_concat_chr${chr}.txt"
    
    bcftools concat -n -f "$OUTPUT_DIR/phase_rare/list_concat_chr${chr}.txt" -Ob -o "$OUTPUT_DIR/ukb_c${chr}_phased.bcf"
    bcftools index "$OUTPUT_DIR/ukb_c${chr}_phased.bcf"
}

# Function to process chromosome X
process_chrX() {
    log "Processing chromosome X"
    
    # Remove PAR regions
    bcftools view -r X:2781480-155701384 "$OUTPUT_DIR/QC/ukb_cX_qc.bcf" -Ob -o "$OUTPUT_DIR/QC/ukb_cX_nonPAR.bcf"
    bcftools index "$OUTPUT_DIR/QC/ukb_cX_nonPAR.bcf"
    
    # Phase common variants
    shapeit5_phase_common --input "$OUTPUT_DIR/QC/ukb_cX_nonPAR.bcf" \
                          --map "$MAP_DIR/chrX.b38.gmap.gz" \
                          --output "$OUTPUT_DIR/phase_common/ukb_cX_common.bcf" \
                          --thread $THREADS \
                          --filter-maf 0.001 \
                          --pedigree "$PED_FILE" \
                          --haploid /path/to/haploid_males.txt
    
    # Ligate common variants
    shapeit5_ligate --input "$OUTPUT_DIR/phase_common/ukb_cX_common.bcf" \
                    --output "$OUTPUT_DIR/ligate/ukb_cX_common_ligated.bcf" \
                    --thread $THREADS
    
    # Phase rare variants
    shapeit5_phase_rare --input "$OUTPUT_DIR/QC/ukb_cX_nonPAR.bcf" \
                        --input-scaffold "$OUTPUT_DIR/ligate/ukb_cX_common_ligated.bcf" \
                        --output "$OUTPUT_DIR/phase_rare/ukb_cX_rare.bcf" \
                        --thread $THREADS \
                        --haploid /path/to/haploid_males.txt
}

# Function to clean up intermediate files
cleanup() {
    log "Cleaning up intermediate files"
    rm -rf "$OUTPUT_DIR/QC" "$OUTPUT_DIR/phase_common" "$OUTPUT_DIR/ligate" "$OUTPUT_DIR/chunks"
}

# Main pipeline
main() {
    log "Starting SHAPEIT5 phasing pipeline"
    
    install_dependencies
    split_vcfs
    
    # Process autosomes
    for chr in {1..22}; do
        while read chunk; do
            perform_qc $chr $chunk
            phase_common $chr $chunk
        done < "$MAP_DIR/chunks_chr${chr}.txt"
        
        ligate_common $chr
        
        while read chunk; do
            phase_rare $chr $chunk
        done < "$MAP_DIR/chunks_chr${chr}.txt"
        
        concatenate_rare $chr
    done
    
    # Process X chromosome
    process_chrX
    
    cleanup
    
    log "Pipeline completed successfully"
}

# Run the main pipeline
main

# TODO: Implement parallelization using GNU Parallel or job arrays for cluster environments
# TODO: Optimize resource management (CPU, memory) based on the specific computational environment
