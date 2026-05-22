#!/bin/bash
# Uso: bash map_reads.sh <SRA_ID> <ARQUIVO_FASTA_REFERENCIA>

SID=$1
REF_FASTA=$2

echo ">>> Iniciando Mapeamento de Reads para $SID..."

# 1. Trava de Segurança: Verifica se a pasta do SRA existe
if [ ! -d "$SID" ]; then
    echo ">>> [ERRO BASH] A pasta do projeto ($SID) não foi encontrada!"
    exit 1
fi

cd "$SID" || exit
mkdir -p 06_Mapping_Results
cd 06_Mapping_Results || exit

# 2. Cria o índice usando o caminho absoluto que o Python mandou (sem o ../../)
echo "Construindo índice do Bowtie2..."
bowtie2-build "$REF_FASTA" alvo_index

# 3. Trava de Segurança: Verifica se as reads existem antes de alinhar
echo "Mapeando reads contra o alvo..."
if [ -f "../${SID}_2.fastq.gz" ]; then
    # Paired-end
    bowtie2 -x alvo_index -1 "../${SID}_1.fastq.gz" -2 "../${SID}_2.fastq.gz" -S alinhamento.sam
elif [ -f "../${SID}.fastq.gz" ]; then
    # Single-end
    bowtie2 -x alvo_index -U "../${SID}.fastq.gz" -S alinhamento.sam
else
    echo ">>> [ERRO BASH] Arquivos brutos .fastq.gz não encontrados na pasta $SID!"
    exit 1
fi

# 4. Processa e conta
echo "Processando arquivos com Samtools..."
samtools view -bS alinhamento.sam > alinhamento.bam
samtools sort alinhamento.bam -o alinhamento_sorted.bam
samtools index alinhamento_sorted.bam

echo "Gerando contagem de reads (idxstats)..."
samtools idxstats alinhamento_sorted.bam > expressao_reads.txt

# Limpeza
rm alinhamento.sam alinhamento.bam
echo ">>> Mapeamento concluído! Resultado salvo em: $SID/06_Mapping_Results/expressao_reads.txt"
