#!/bin/bash

# A bash script to automate *de novo* transcriptome assembly using Trinity.
# The script will stop immediately if any command fails.
set -e

## --- Configuration ---
# Path to the Trimmomatic adapter file.
# You may need to change this to the correct path on your system.
ADAPTERS="/usr/share/trimmomatic/TruSeq3-PE-2.fa"
# Maximum RAM (in Gb) to use for Trinity assembly.
# Set to 8 for low-resource systems, adjust if needed (e.g., 4 or 16).
MAX_RAM=8
# Number of CPU threads to use
THREADS=2

## --- Input Validation ---
if [ "$#" -ne 2 ]; then
    echo "ERROR: You must provide two arguments."
    echo "Usage: $0 <forward_reads.fastq.gz> <reverse_reads.fastq.gz>"
    echo "Example: $0 SRR31116834_1.fastq.gz SRR31116834_2.fastq.gz"
    exit 1
fi

# Assign input files to variables for clarity
READ1=$1
READ2=$2

# Create a base name for output files
PREFIX=$(basename ${READ1} _1.fastq.gz)
echo "=> Starting Trinity pipeline for sample: ${PREFIX}"

## --- Step 1: Raw Read Quality Control (FastQC) ---
echo "### STEP 1: Running FastQC on raw reads ###"
mkdir -p 01_QC_Reports
fastqc --threads ${THREADS} -o 01_QC_Reports ${READ1} ${READ2}

## --- Step 2: Adapter and Quality Trimming (Trimmomatic) ---
echo "### STEP 2: Trimming adapters and low-quality reads with Trimmomatic ###"
mkdir -p 02_Trimmed_Reads
TrimmomaticPE -threads ${THREADS} -phred33 \
    ${READ1} ${READ2} \
    02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_1_unpaired.fastq.gz \
    02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

## --- Step 3: De Novo Assembly (Trinity) ---
# We do not use bbnorm here because Trinity has built-in in silico read normalization, which it runs by default
# to keep memory usage low.
echo "### STEP 3: Assembling trimmed paired reads with Trinity ###"
OUT_DIR="03_Trinity_Assembly_${PREFIX}"

Trinity --seqType fq \
    --left 02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz \
    --right 02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz \
    --CPU ${THREADS} \
    --max_memory ${MAX_RAM}G \
    --output ${OUT_DIR}

## --- Step 4: Post-Processing QC and Aggregation (FastQC + MultiQC) ---
echo "### STEP 4: Running FastQC on trimmed reads and creating MultiQC report ###"
# Run FastQC on the trimmed paired reads
fastqc --threads ${THREADS} -o 01_QC_Reports 02_Trimmed_Reads/${PREFIX}_*_paired.fastq.gz

# Run MultiQC to aggregate all reports
echo "=> Aggregating QC reports with MultiQC..."
mkdir -p 04_MultiQC_Report
multiqc . -o 04_MultiQC_Report

## --- Pipeline Complete ---
echo "✅ Pipeline finished successfully!"
echo "=> Final QC report is in '04_MultiQC_Report/multiqc_report.html'"
echo "=> Assembly results are in the '${OUT_DIR}' directory."
echo "=> The final assembled sequences are in '${OUT_DIR}/Trinity.fasta'"
