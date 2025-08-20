#!/bin/bash
# Process mapped reads for a trio: father, mother, child
# Usage: bash process_trio.sh

# Paths to SAM files
FATHER_SAM=~/OneDrive/Om/exomesequencing/mapped/father.sam
MOTHER_SAM=~/OneDrive/Om/exomesequencing/mapped/mother.sam
CHILD_SAM=~/OneDrive/Om/exomesequencing/mapped/proband.sam

# Output directory
OUTDIR=~/OneDrive/Om/exomesequencing/mapped/processed
mkdir -p $OUTDIR

# Function to process a sample
process_sample() {
    SAMPLE_NAME=$1
    SAM_FILE=$2

    echo "Processing $SAMPLE_NAME ..."

    # 1️⃣ Convert SAM → BAM and sort by coordinate
    samtools view -bS $SAM_FILE > $OUTDIR/${SAMPLE_NAME}.bam

    # 2️⃣ Filter for properly paired reads (both mates mapped)
    samtools view -b -f 3 $OUTDIR/${SAMPLE_NAME}.bam > $OUTDIR/${SAMPLE_NAME}.filtered.bam

    # 3️⃣ Sort by queryname for fixmate
    samtools sort -n -o $OUTDIR/${SAMPLE_NAME}.qname.sorted.bam $OUTDIR/${SAMPLE_NAME}.filtered.bam

    # 4️⃣ Run fixmate
    samtools fixmate -m $OUTDIR/${SAMPLE_NAME}.qname.sorted.bam $OUTDIR/${SAMPLE_NAME}.fixmate.bam

    # 5️⃣ Sort by coordinate for markdup
    samtools sort -o $OUTDIR/${SAMPLE_NAME}.coord.sorted.bam $OUTDIR/${SAMPLE_NAME}.fixmate.bam

    # 6️⃣ Remove duplicates
    samtools markdup -r $OUTDIR/${SAMPLE_NAME}.coord.sorted.bam $OUTDIR/${SAMPLE_NAME}.dedup.bam

    # 7️⃣ Index final BAM
    samtools index $OUTDIR/${SAMPLE_NAME}.dedup.bam

    echo "$SAMPLE_NAME done."
}

# Process all three samples
process_sample father $FATHER_SAM
process_sample mother $MOTHER_SAM
process_sample child $CHILD_SAM

echo "All samples processed successfully!"
