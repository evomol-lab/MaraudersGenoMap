#!/bin/bash

# A bash script to set up the environment for the genome assembly pipeline.
# The script will stop immediately if any command fails.
set -e

## --- 0. Environment Setup ---
echo "Installing necessary tools if not already installed..."
sudo apt-get update
sudo apt-get install -y sra-toolkit pigz coreutils fastqc parallel
sudo apt-get install -y fastqc trimmomatic bbmap multiqc prodigal hmmer seqtk Trinity megahit spades
sudo apt-get install -y python3-pip samtools bamtools
wget http://ftp.br.debian.org/debian/pool/main/n/ncurses/libtinfo5_6.2+20201114-2+deb11u2_amd64.deb
wget http://ftp.br.debian.org/debian/pool/main/n/ncurses/libncurses5_6.2+20201114-2+deb11u2_amd64.deb
sudo dpkg -i libtinfo5_6.2*.deb libncurses5_6.2*.deb
rm libtinfo5_6.2+20201114-2+deb11u2_amd64.deb
rm libncurses5_6.2+20201114-2+deb11u2_amd64.deb


# Choose one of the following methods to install Python packages:
# Method 1: Using pip3 (system-wide installation)
# Method 2: Using pipx (isolated environment installation)
# Uncomment the method you prefer.

#pip3 install biopython multitqc
pip install biopython multiqc
#pipx install biopython multitqc

echo "✅ All required tools are installed."
