#!/bin/bash

# A bash script to set up the environment for the genome assembly pipeline.
# The script will stop immediately if any command fails.
set -e

## --- 0. Environment Setup ---
echo "Installing necessary tools if not already installed..."
sudo apt-get update
sudo apt-get install -y sra-toolkit pigz coreutils md5sum fastqc gnu-parallel
sudo apt-get install -y fastqc trimmomatic bbmap multiqc prodigal hmmer seqtk trinity megahit spades
sudo apt-get install -y python3-pip samtools bamtools

# Choose one of the following methods to install Python packages:
# Method 1: Using pip3 (system-wide installation)
# Method 2: Using pipx (isolated environment installation)
# Uncomment the method you prefer.

#pip3 install biopython multitqc
pip install biopython multitqc
#pipx install biopython multitqc

echo "✅ All required tools are installed."