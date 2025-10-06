#!/bin/bash

# A bash script to automate a *de novo* genome assembly pipeline.
# The script will stop immediately if any command fails.
set -e

## --- Configuration ---
# Path to the Trimmomatic adapter file.
# You may need to change this to the correct path on your system.
ADAPTERS="/usr/share/trimmomatic/TruSeq3-PE-2.fa"

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
echo "=> Starting pipeline for sample: ${PREFIX}"

## --- Step 1: Raw Read Quality Control (FastQC) ---
echo "### STEP 1: Running FastQC on raw reads ###"
mkdir -p 01_QC_Reports
fastqc --threads 2 -o 01_QC_Reports ${READ1} ${READ2}

## --- Step 2: Adapter and Quality Trimming (Trimmomatic) ---
echo "### STEP 2: Trimming adapters and low-quality reads with Trimmomatic ###"
mkdir -p 02_Trimmed_Reads
TrimmomaticPE -phred33 \
    ${READ1} ${READ2} \
    02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_1_unpaired.fastq.gz \
    02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

## --- Step 3: Digital Normalization (bbnorm.sh) ---
echo "### STEP 3: Normalizing read coverage with bbnorm.sh ###"
mkdir -p 03_Normalized_Reads
bbnorm.sh in1=02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz \
          in2=02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz \
          out1=03_Normalized_Reads/${PREFIX}_1_normalized.fastq.gz \
          out2=03_Normalized_Reads/${PREFIX}_2_normalized.fastq.gz \
          target=100

## --- Step 4: Post-Processing QC and Aggregation (FastQC + MultiQC) ---
echo "### STEP 4: Running FastQC on normalized reads and creating MultiQC report ###"
# Run FastQC on the newly normalized paired reads
fastqc -o 01_QC_Reports 03_Normalized_Reads/${PREFIX}_*_normalized.fastq.gz

# Run MultiQC to aggregate all reports
echo "=> Aggregating QC reports with MultiQC..."
mkdir -p 04_MultiQC_Report
multiqc . -o 04_MultiQC_Report

## --- Step 5: De Novo Assembly (SPAdes) ---
echo "### STEP 5: Assembling normalized reads with SPAdes ###"
spades.py --careful --phred-offset 33 --memory 30 \
    -1 03_Normalized_Reads/${PREFIX}_1_normalized.fastq.gz \
    -2 03_Normalized_Reads/${PREFIX}_2_normalized.fastq.gz \
    -o 05_SPAdes_Assembly_${PREFIX}

## --- Pipeline Complete ---
echo "✅ Pipeline finished successfully!"
echo "=> Final QC report is in '04_MultiQC_Report/multiqc_report.html'"
echo "=> Assembly results are in the '05_SPAdes_Assembly_${PREFIX}' directory."