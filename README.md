# **Exome Sequencing Analysis Pipeline**
## Overview

This guide covers downloading sequencing data, preparing the reference genome, aligning reads with `bwa mem`, and viewing the resulting BAM files using `samtools`.

---

## **1️⃣ Download FASTQ Files**

Run this command to download all paired-end FASTQ files:

```bash
wget https://zenodo.org/record/3243160/files/father_R1.fq.gz \
     https://zenodo.org/record/3243160/files/father_R2.fq.gz \
     https://zenodo.org/record/3243160/files/mother_R1.fq.gz \
     https://zenodo.org/record/3243160/files/mother_R2.fq.gz \
     https://zenodo.org/record/3243160/files/proband_R1.fq.gz \
     https://zenodo.org/record/3243160/files/proband_R2.fq.gz
```

---

# Create the Conda environment from the YAML file
```
conda env create -f environment.yml
```

# Activate the environment
```
conda activate exome-seq
```

## **3️⃣ Prepare Reference Genome**

1. Decompress if needed:

```bash
gunzip hg19_chr8.fa.gz
```

2. Index for BWA:

```bash
bwa index hg19_chr8.fa
```

3. Index for SAMtools:

```bash
samtools faidx hg19_chr8.fa
```

---

## **4️⃣ Align Reads with BWA MEM**

```bash
# Father
bwa mem -R '@RG\tID:000\tSM:father\tPL:ILLUMINA' \
hg19_chr8.fa father_R1.fq.gz father_R2.fq.gz > father.sam

# Mother
bwa mem -R '@RG\tID:001\tSM:mother\tPL:ILLUMINA' \
hg19_chr8.fa mother_R1.fq.gz mother_R2.fq.gz > mother.sam

# Proband
bwa mem -R '@RG\tID:002\tSM:proband\tPL:ILLUMINA' \
hg19_chr8.fa proband_R1.fq.gz proband_R2.fq.gz > proband.sam
```

---

## **5️⃣ Process SAM → BAM → Deduplicated BAMs**

Here’s a **one-shot script** to convert, filter, fixmate, sort, mark duplicates, and index:

```bash
#!/bin/bash
OUTDIR=mapped/processed
mkdir -p $OUTDIR

process_sample() {
    SAMPLE=$1
    SAM=$2
    echo "Processing $SAMPLE ..."

    samtools view -bS $SAM > $OUTDIR/${SAMPLE}.bam
    samtools view -b -f 3 $OUTDIR/${SAMPLE}.bam > $OUTDIR/${SAMPLE}.filtered.bam
    samtools sort -n -o $OUTDIR/${SAMPLE}.qname.sorted.bam $OUTDIR/${SAMPLE}.filtered.bam
    samtools fixmate -m $OUTDIR/${SAMPLE}.qname.sorted.bam $OUTDIR/${SAMPLE}.fixmate.bam
    samtools sort -o $OUTDIR/${SAMPLE}.coord.sorted.bam $OUTDIR/${SAMPLE}.fixmate.bam
    samtools markdup -r $OUTDIR/${SAMPLE}.coord.sorted.bam $OUTDIR/${SAMPLE}.dedup.bam
    samtools index $OUTDIR/${SAMPLE}.dedup.bam
    echo "$SAMPLE done."
}

process_sample father father.sam
process_sample mother mother.sam
process_sample proband proband.sam
```

*Save as `process_trio.sh`, make executable:*

```bash
chmod +x process_trio.sh
./process_trio.sh
```

✅ Output:

* `father.dedup.bam`, `mother.dedup.bam`, `proband.dedup.bam` (indexed)

---

## **6️⃣ Variant Calling with FreeBayes**

```bash
freebayes -f hg19_chr8.fa father.dedup.bam mother.dedup.bam proband.dedup.bam > trio_variants.vcf
```

*Optional parameters for simple diploid calling:*

```bash
freebayes -f hg19_chr8.fa --ploidy 2 --min-alternate-count 1 --min-alternate-fraction 0.2 \
father.dedup.bam mother.dedup.bam proband.dedup.bam > trio_variants.vcf
```

---

## **7️⃣ Normalize VCF (Left-align & Split multiallelics)**

```bash
# Prepend contig header
echo "##contig=<ID=chr8,length=$(cut -f2 hg19_chr8.fa.fai)>" > contig_header.txt
bcftools reheader -h <(cat contig_header.txt <(grep '^#' trio_variants.vcf)) \
-o trio_variants.headered.vcf trio_variants.vcf

# Normalize
bcftools norm -f hg19_chr8.fa -m -both --check-ref w \
-o trio_variants.norm.vcf trio_variants.headered.vcf

# Optional compression & indexing
bgzip trio_variants.norm.vcf
tabix -p vcf trio_variants.norm.vcf.gz
```

---

## **8️⃣ Annotate with SnpEff**

```bash
# Download GRCh37.75 database
java -jar snpEff.jar download GRCh37.75

# Annotate VCF
java -jar snpEff.jar -v GRCh37.75 trio_variants.norm.vcf > trio_variants.snpeff.vcf
```

---

## **9️⃣ Load Annotated VCF into GEMINI**

```bash
gemini load -v trio_variants.snpeff.vcf -p trio.ped --cores 4 --skip-sanity-check --force trio.db
```

*Check variants:*

```bash
gemini query -q "SELECT count(*) FROM variants;" trio.db
gemini de_novo trio.db
```

---

## **10️⃣ Optional: Inspect & Export**

```bash
# Filter by impact
gemini query -q "SELECT * FROM variants WHERE impact_severity='HIGH' OR impact_severity='MED';" trio.db

# Export CSV
gemini query -q "SELECT chrom, start, end, gene, impact, genotype_father, genotype_mother, genotype_proband FROM variants;" trio.db -o annotated_variants.csv
```

---

✅ **Outcome:**

* Deduplicated, indexed BAMs
* Normalized VCF
* SnpEff-annotated VCF
* GEMINI database ready for queries, inheritance analysis, and reporting

---

