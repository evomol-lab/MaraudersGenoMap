#!/bin/bash

# A script to extract protein sequences from a FASTA file based on a HMMER domain table output.
# It will stop immediately if any command fails.
set -e

## --- Input Validation ---
if [ "$#" -ne 2 ]; then
    echo "ERROR: You must provide two arguments."
    echo "Usage:   $0 <hmmer_table.tbl> <proteins.faa>"
    echo "Example: $0 PF00139_domain_hits.tbl predicted_proteins.faa"
    exit 1
fi

# Assign input files to variables
HMMER_TBL=$1
PROTEINS_FAA=$2

# Define output filenames
ID_LIST="protein_ids_to_extract.txt"
OUTPUT_FASTA="lectin_hits.faa"

## --- Step 1: Create a Unique List of Protein IDs ---
echo "### STEP 1: Parsing HMMER output to get unique protein IDs... ###"

# 1. `grep -v '^#'`: Removes all comment lines from the HMMER output.
# 2. `awk '{print $1}'`: Prints only the first column, which is the protein ID.
# 3. `sort -u`: Sorts the IDs and removes duplicates.
grep -v '^#' "${HMMER_TBL}" | awk '{print $1}' | sort -u > "${ID_LIST}"

# Count how many unique proteins were found
COUNT=$(wc -l < "${ID_LIST}")
echo "=> Found ${COUNT} unique proteins with lectin domains."

## --- Step 2: Extract FASTA Sequences with seqtk ---
echo "### STEP 2: Extracting full protein sequences using seqtk... ###"
seqtk subseq "${PROTEINS_FAA}" "${ID_LIST}" > "${OUTPUT_FASTA}"

## --- Cleanup and Final Report ---
# Remove the intermediate list of IDs
rm "${ID_LIST}"

echo "✅ Success!"
echo "=> The ${COUNT} full protein sequences have been saved to: ${OUTPUT_FASTA}"