#!/bin/bash
# Shell script to download and compress SRA data for paired-end reads from Leguminosae species
# Usage: ./get.sh <SRA_accession>
# Example: ./get.sh SRR1234567 # Replace SRR1234567 with the actual SRA accession number
# Note: Ensure you have 'prefetch' and 'fasterq-dump' from the SRA Toolkit installed and configured.
# Note: Ensure you have 'pigz' installed for parallel gzip compression.
# Note: Ensure you have 'md5sum' installed for checksum verification.

# It stops immediately if any command fails.
set -e

# --- 0. Environment Setup ---
# Uncomment the following lines if you need to install the required tools.
# You may need sudo privileges to install these tools.
#echo "Installing necessary tools if not already installed..."
#sudo apt-get update
#sudo apt-get install -y sra-toolkit pigz coreutils md5sum

# --- 1. Input Validation ---
# Check if an SRA accession number was provided as an argument.
if [ -z "$1" ]; then
  echo "Usage: $0 <SRA_ACCESSION>"
  echo "Example: $0 SRR31116834"
  exit 1
fi

SRA_ACCESSION=$1
echo "Processing SRA Accession: ${SRA_ACCESSION}"

# --- 2. Download SRA Data ---
echo "=> Step 1: Prefetching SRA data..."
prefetch --progress "${SRA_ACCESSION}"

# --- 3. Extract FASTQ files ---
echo "=> Step 2: Extracting FASTQ files using fasterq-dump..."
fasterq-dump --split-files --progress "${SRA_ACCESSION}"

# --- 4. Compress FASTQ files with pigz ---
echo "=> Step 3: Compressing FASTQ files with pigz..."
# pigz will automatically use multiple cores for faster compression.
# It also removes the original file by default upon successful compression.
if [ -f "${SRA_ACCESSION}_1.fastq" ] && [ -f "${SRA_ACCESSION}_2.fastq" ]; then
  pigz "${SRA_ACCESSION}_1.fastq"
  pigz "${SRA_ACCESSION}_2.fastq"
  echo "=> Compression successful."
else
  echo "Error: FASTQ files not found after dump. Cannot compress."
  exit 1
fi

# --- 5. Final Report ---
echo "✅ Done!"
echo "Output files:"
echo "- ${SRA_ACCESSION}_1.fastq.gz"
echo "- ${SRA_ACCESSION}_2.fastq.gz"