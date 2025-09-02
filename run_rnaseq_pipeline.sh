#!/bin/bash

# ----------------------------
# Automated RNA-Seq Pipeline for HB3 Plasmodium (Combined FASTQs)
# Tools: HISAT2, SAMtools, featureCounts
# ----------------------------

# === CONFIG ===
GENOME_FASTA="PlasmoDB-61_PfalciparumHB3_Genome.fasta"
ANNOTATION_GFF="PlasmoDB-61_PfalciparumHB3.gff"
INDEX_DIR="Genomes_HB3"
INDEX_NAME="HB3indexed"
THREADS=12
OUTDIR="pipeline_output"

# === PREP ===
echo "[INFO] Checking required tools..."
for tool in hisat2 samtools featureCounts; do
    if ! command -v $tool &> /dev/null; then
        echo "[ERROR] $tool not found in PATH. Please install it first."
        exit 1
    fi
done

mkdir -p $OUTDIR/{bam,logs,counts}

# === BUILD INDEX ===
if [ ! -f ${INDEX_DIR}/${INDEX_NAME}.1.ht2 ]; then
    echo "[INFO] Building HISAT2 index..."
    mkdir -p $INDEX_DIR
    hisat2-build $GENOME_FASTA ${INDEX_DIR}/${INDEX_NAME} | tee $OUTDIR/logs/index_build.log
else
    echo "[INFO] HISAT2 index already exists. Skipping."
fi

# === COMBINE FASTQ FILES ===
echo "[INFO] Combining R1 files..."
cat N175-32_S32_L001_R1_001.fastq.gz N175-32_S32_L002_R1_001.fastq.gz N175-32_S32_L003_R1_001.fastq.gz N175-32_S32_L004_R1_001.fastq.gz > combined_R1.fastq.gz

echo "[INFO] Combining R2 files..."
cat N175-32_S32_L001_R2_001.fastq.gz N175-32_S32_L002_R2_001.fastq.gz N175-32_S32_L003_R2_001.fastq.gz N175-32_S32_L004_R2_001.fastq.gz > combined_R2.fastq.gz

# === ALIGN COMBINED SAMPLE ===
echo "[INFO] Aligning combined sample..."
hisat2 -p 12 --max-intronlen 10000 -x ${INDEX_DIR}/${INDEX_NAME} \
    -1 combined_R1.fastq.gz -2 combined_R2.fastq.gz 2> $OUTDIR/logs/combined_hisat2.log | \
    samtools view -Sb - | \
    samtools sort -o $OUTDIR/bam/combined.bam

samtools flagstat $OUTDIR/bam/combined.bam > $OUTDIR/logs/combined_mapping_stats.txt
samtools index $OUTDIR/bam/combined.bam

# === RUN FEATURECOUNTS ===
echo "[INFO] Running featureCounts..."
featureCounts -T $THREADS -p -O -Q 5 -t exon -g ID \
    -a $ANNOTATION_GFF \
    -o $OUTDIR/counts/combined_counts.txt $OUTDIR/bam/combined.bam \
    2> $OUTDIR/logs/featureCounts.log

echo "[INFO] Pipeline completed! All outputs are in the $OUTDIR directory."

