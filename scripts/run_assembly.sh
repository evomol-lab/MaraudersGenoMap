#!/bin/bash
set -e

# Argumentos: R1, R2, Threads, RAM, Assembler, DelRaw, DelTrim
READ1=$1; READ2=$2; THREADS=$3; MEMORY=$4; ASSEMBLER=$5; DEL_RAW=$6; DEL_TRIM=$7
MEMORY_BYTES=$(( MEMORY * 1024 * 1024 * 1024 ))
READ1_NAME=$(basename "$READ1")
if [[ "$READ1_NAME" == *_1.fastq.gz ]]; then
    PREFIX="${READ1_NAME%_1.fastq.gz}"
elif [[ "$READ1_NAME" == *.fastq.gz ]]; then
    PREFIX="${READ1_NAME%.fastq.gz}"
else
    PREFIX="${READ1_NAME%.*}"
fi
ADAPTERS="/usr/share/trimmomatic/TruSeq3-PE-2.fa"

mkdir -p 01_QC_Reports 02_Trimmed_Reads 03_Normalized_Reads 04_MultiQC_Report 05_Assembly_Results

echo "### [STEP 1] 📊 QC Inicial ###"
if [ "$READ2" != "none" ] && [ -f "$READ2" ]; then
    fastqc --threads 2 -o 01_QC_Reports "$READ1" "$READ2"
else
    fastqc --threads 2 -o 01_QC_Reports "$READ1"
fi

echo "### [STEP 2] ✂️  Trimming ###"
if [ "$READ2" != "none" ] && [ -f "$READ2" ]; then
    TrimmomaticPE -threads "$THREADS" -phred33 "$READ1" "$READ2" \
        02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_1_unpaired.fastq.gz \
        02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz 02_Trimmed_Reads/${PREFIX}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    R1_FINAL="02_Trimmed_Reads/${PREFIX}_1_paired.fastq.gz"
    R2_FINAL="02_Trimmed_Reads/${PREFIX}_2_paired.fastq.gz"
else
    TrimmomaticSE -threads "$THREADS" -phred33 "$READ1" 02_Trimmed_Reads/${PREFIX}_trimmed.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    R1_FINAL="02_Trimmed_Reads/${PREFIX}_trimmed.fastq.gz"
fi

echo "### [STEP 3] 🧬 Normalização ###"
if [[ "$ASSEMBLER" != "trinity" ]]; then
    if [ "$READ2" != "none" ] && [ -f "$READ2" ]; then
        bbnorm.sh in1="$R1_FINAL" in2="$R2_FINAL" out1=03_Normalized_Reads/${PREFIX}_1_norm.fastq.gz out2=03_Normalized_Reads/${PREFIX}_2_norm.fastq.gz target=100
        R1_ASM="03_Normalized_Reads/${PREFIX}_1_norm.fastq.gz"
        R2_ASM="03_Normalized_Reads/${PREFIX}_2_norm.fastq.gz"
    else
        bbnorm.sh in="$R1_FINAL" out=03_Normalized_Reads/${PREFIX}_norm.fastq.gz target=100
        R1_ASM="03_Normalized_Reads/${PREFIX}_norm.fastq.gz"
    fi
else
    R1_ASM="$R1_FINAL"; R2_ASM="$R2_FINAL"
fi

echo "### [STEP 4] 📈 MultiQC ###"
multiqc . -o 04_MultiQC_Report

echo "### [STEP 5] 🏗️  Assembly Final com $ASSEMBLER ###"
case $ASSEMBLER in
    megahit)
        megahit $( [ "$READ2" != "none" ] && echo "-1 $R1_ASM -2 $R2_ASM" || echo "-r $R1_ASM" ) -t "$THREADS" -m "$MEMORY_BYTES" -o "05_Assembly_Results/MEGAHIT_${PREFIX}" ;;
    spades)
        spades.py --careful -t "$THREADS" -m "$MEMORY" -o "05_Assembly_Results/SPADES_${PREFIX}" $( [ "$READ2" != "none" ] && echo "-1 $R1_ASM -2 $R2_ASM" || echo "-s $R1_ASM" ) ;;
    trinity)
        Trinity --seqType fq --max_memory "${MEMORY}G" --CPU "$THREADS" --output "05_Assembly_Results/TRINITY_${PREFIX}" $( [ "$READ2" != "none" ] && echo "--left $R1_ASM --right $R2_ASM" || echo "--single $R1_ASM" ) ;;
esac

# Limpeza automática baseada na seleção do Dashboard (Sem perguntas!)
[ "$DEL_RAW" == "s" ] && rm -f "${PREFIX}"_*.fastq.gz
[ "$DEL_TRIM" == "s" ] && rm -rf 02_Trimmed_Reads

echo -e "\a\n✅ Pipeline Marauders finalizado com sucesso!"
