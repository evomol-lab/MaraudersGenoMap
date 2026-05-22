#!/bin/bash
set -e

SRA_ACCESSION=$1
LAYOUT=$2 # 'single' ou 'paired'

echo -e "\n=> Iniciando download para ID: $SRA_ACCESSION | Layout: $LAYOUT"

# --- ORGANIZAÇÃO CENTRALIZADA ---
mkdir -p "$SRA_ACCESSION"
cd "$SRA_ACCESSION"

# Reset de segurança para evitar travamentos de rede
vdb-config --restore > /dev/null 2>&1 || true
rm -f *.lock 2>/dev/null || true

# Download (Lógica para não travar no SRA Lite)
echo "=> Passo 1/3: Baixando dados (prefetch)..."
if [ "$LAYOUT" == "single" ]; then
    prefetch --progress --type sra "${SRA_ACCESSION}"
else
    prefetch --progress "${SRA_ACCESSION}"
fi

# Extração
echo "=> Passo 2/3: Extraindo FASTQ..."
if [ "$LAYOUT" == "single" ]; then
    if ! fasterq-dump --split-files --progress "${SRA_ACCESSION}"; then
        echo "=> Fallback para fastq-dump..."
        fastq-dump --skip-technical "${SRA_ACCESSION}"
    fi
else
    fasterq-dump --split-3 --progress "${SRA_ACCESSION}"
fi

# Compactação
echo "=> Passo 3/3: Compactando com pigz..."
sleep 2
pigz -f ${SRA_ACCESSION}*.fastq

echo "✅ Ambiente pronto na pasta: $SRA_ACCESSION"
