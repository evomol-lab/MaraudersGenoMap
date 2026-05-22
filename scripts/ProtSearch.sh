#!/bin/bash
set -e

## --- Validação de Entrada ---
# Agora aceitamos o 3º argumento para o MODO (single ou meta)
if [ "$#" -ne 3 ]; then
    echo "ERRO: Você deve fornecer três argumentos."
    echo "Uso:   $0 <contigs.fasta> <profile.hmm> <modo: single|meta>"
    exit 1
fi

CONTIGS=$1
HMM_PROFILE=$2
MODE=$3

# Nome base para os resultados
PREFIX=$(basename ${HMM_PROFILE} .hmm)
echo "=> Iniciando busca de proteínas para a família: ${PREFIX}"
echo "=> Modo Prodigal: ${MODE}"

## --- Passo 1: Predizer Proteínas com Prodigal ---
echo "### PASSO 1: Predizendo genes com Prodigal ###"
mkdir -p 01_Predicted_Proteins
PROTEINS="01_Predicted_Proteins/predicted_proteins.faa"

# -q: modo silencioso (evita travar o Dashboard)
# -o /dev/null: descarta a saída de texto gigante no terminal
prodigal -i "${CONTIGS}" -a "${PROTEINS}" -p "${MODE}" -q -o /dev/null

## --- Passo 2: Busca com HMMER ---
echo "### PASSO 2: Buscando domínios com HMMER ###"
mkdir -p 02_HMMER_Results
HMM_OUTPUT="02_HMMER_Results/${PREFIX}_domain_hits.tbl"

# Executa a busca
hmmsearch --domtblout "${HMM_OUTPUT}" "${HMM_PROFILE}" "${PROTEINS}"

echo "✅ Busca concluída!"
echo "=> Proteínas em: ${PROTEINS}"
echo "=> Tabela de hits em: ${HMM_OUTPUT}"
