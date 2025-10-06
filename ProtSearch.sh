#!/bin/bash

# A script to predict proteins from assembled contigs and search for a protein family using HMMER.
# The script will stop immediately if any command fails.
set -e

## --- Input Validation ---
if [ "$#" -ne 2 ]; then
    echo "ERROR: You must provide two arguments."
    echo "Usage:   $0 <contigs.fasta> <profile.hmm>"
    echo "Example: $0 05_MEGAHIT_Assembly/final.contigs.fa PF00139.hmm"
    exit 1
fi

# Assign input files to variables
CONTIGS=$1
HMM_PROFILE=$2

# Create a base name for the output file from the HMM profile name
PREFIX=$(basename ${HMM_PROFILE} .hmm)
echo "=> Starting protein search for family: ${PREFIX}"

## --- Step 1: Predict Proteins with Prodigal ---
echo "### STEP 1: Predicting protein-coding genes with Prodigal ###"
mkdir -p 01_Predicted_Proteins
PROTEINS="01_Predicted_Proteins/predicted_proteins.faa"
prodigal -i ${CONTIGS} -a ${PROTEINS} -p meta

## --- Step 2: Search for Protein Family with HMMER ---
echo "### STEP 2: Searching for ${PREFIX} domains with HMMER ###"
mkdir -p 02_HMMER_Results
HMM_OUTPUT="02_HMMER_Results/${PREFIX}_domain_hits.tbl"
hmmsearch --domtblout ${HMM_OUTPUT} ${HMM_PROFILE} ${PROTEINS}

## --- Pipeline Complete ---
echo "✅ Search complete!"
echo "=> Predicted proteins are in: ${PROTEINS}"
echo "=> HMMER domain search results are in: ${HMM_OUTPUT}"