
# Exome Sequencing Analysis Pipeline

## Overview

This guide covers downloading sequencing data, preparing the reference genome, aligning reads with `bwa mem`, and viewing the resulting BAM files using `samtools`.

---

## 1. Download FASTQ Files

Run this command to download all paired-end FASTQ files:

```
wget https://zenodo.org/record/3243160/files/father_R1.fq.gz \
     https://zenodo.org/record/3243160/files/father_R2.fq.gz \
     https://zenodo.org/record/3243160/files/mother_R1.fq.gz \
     https://zenodo.org/record/3243160/files/mother_R2.fq.gz \
     https://zenodo.org/record/3243160/files/proband_R1.fq.gz \
     https://zenodo.org/record/3243160/files/proband_R2.fq.gz
````

---

## 2. Install Required Software

### Using Conda

```
conda install -c bioconda bwa samtools
```

### Manual Installation of BWA (optional)

```
git clone https://github.com/lh3/bwa.git
cd bwa
make
sudo cp bwa /usr/local/bin/
```

---

## 3. Prepare Reference Genome

If your reference genome is compressed (`.fa.gz`), decompress it first:

```
gunzip hg19_chr8.fa.gz
```

Index the reference genome for `bwa`:

```
bwa index hg19_chr8.fa
```

---

## 4. Align Reads with BWA MEM

For paired-end reads (recommended):

```
bwa mem hg19_chr8.fa father_R1.fq.gz father_R2.fq.gz > father.sam
bwa mem hg19_chr8.fa mother_R1.fq.gz mother_R2.fq.gz > mother.sam
bwa mem hg19_chr8.fa proband_R1.fq.gz proband_R2.fq.gz > proband.sam
```

For single-end reads (if only one FASTQ file is available):

```
bwa mem hg19_chr8.fa father_R1.fq.gz > father_single.sam
```

---

## 5. Convert SAM to BAM

```
samtools view -Sb father.sam > father.bam
```

Repeat for other samples as needed.

---

## 6. View BAM Files

* View BAM as text:

  ```
  samtools view father.bam | less -S
  ```

* View BAM with header:

  ```
  samtools view -h father.bam | less -S
  ```

* Get alignment statistics:

  ```
  samtools flagstat father.bam
  ```

* Count reads:

  ```
  samtools view -c father.bam
  ```

---

## 7. Exit `less` Viewer

When viewing BAM or SAM output with `less`, press:

```
q
```

to quit and return to the command prompt.

---

## Notes

* Always decompress the reference genome before indexing and alignment.
* `bwa mem` handles paired-end reads by providing both R1 and R2 FASTQ files.
* Use `samtools` for BAM file processing and viewing.

---
## SeqKit Stats
```
seqkit stats *gz
seqkit -less 20 father_R1.fq.gz 
seqkit stats *gz -a
```

## BWA MEM2
```
bwa index hg19.fa

```
### Map with paired ENd for father
```
bwa mem -R '@RG\tID:000\tSM:father\tPL:ILLUMINA' \
~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
~/OneDrive/Om/exomesequencing/raw_data/father_R1.fq.gz \
~/OneDrive/Om/exomesequencing/raw_data/father_R2.fq.gz \
> ~/OneDrive/Om/exomesequencing/mapped/father.sam


```
### Map with paired ENd for mother
```
bwa mem -R '@RG\tID:001\tSM:mother\tPL:ILLUMINA' \
~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
~/OneDrive/Om/exomesequencing/raw_data/mother_R1.fq.gz \
~/OneDrive/Om/exomesequencing/raw_data/mother_R2.fq.gz \
> ~/OneDrive/Om/exomesequencing/mapped/mother.sam

```
### Map with paired ENd for proband
```
bwa mem -R '@RG\tID:002\tSM:proband\tPL:ILLUMINA' \
~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
~/OneDrive/Om/exomesequencing/raw_data/proband_R1.fq.gz \
~/OneDrive/Om/exomesequencing/raw_data/proband_R2.fq.gz \
> ~/OneDrive/Om/exomesequencing/mapped/proband.sam

```
### Father
```
samtools view -bS father.sam | samtools sort -o father.sorted.bam
```

### Mother
```
samtools view -bS mother.sam | samtools sort -o mother.sorted.bam
```

### Child/Proband
```
samtools view -bS proband.sam | samtools sort -o child.sorted.bam
```

### Output
-bS → convert SAM to BAM

## Filter 

### Father
```
samtools sort -n -o father.qname.sorted.bam father.filtered.bam
```

### Mother
```
samtools sort -n -o mother.qname.sorted.bam mother.filtered.bam
```

### Child
```
samtools sort -n -o child.qname.sorted.bam child.filtered.bam
```

### Note
-f 3 → keep only properly paired reads.
# Fixmate
### Father
```
samtools fixmate -m father.filtered.bam father.fixmate.bam
```

### Mother
```
samtools fixmate -m mother.filtered.bam mother.fixmate.bam
```

### Child
```
samtools fixmate -m child.filtered.bam child.fixmate.bam
```

### Short Reads
```
samtools sort -o father.coord.sorted.bam father.fixmate.bam
samtools sort -o mother.coord.sorted.bam mother.fixmate.bam
samtools sort -o child.coord.sorted.bam child.fixmate.bam
```


### 3. Deduplicate reads

```
samtools markdup -r father.coord.sorted.bam father.dedup.bam
samtools markdup -r mother.coord.sorted.bam mother.dedup.bam
samtools markdup -r child.coord.sorted.bam child.dedup.bam

```

### Note: samtools markdup requires coordinate-sorted BAMs. Our filtered.bam files are already sorted from step 1, so this works.

### 4. Index the final BAMs
```
samtools index father.dedup.bam
samtools index mother.dedup.bam
samtools index child.dedup.bam
```
✅ Key points to remember:

fixmate → BAM must be queryname sorted.

markdup → BAM must be coordinate sorted.

Always filter proper pairs (-f 3) before fixmate/markdup.

If you want, I can write a full one-shot  script that does:

SAM → BAM → filter proper pairs

Queryname sort → fixmate

Coordinate sort → markdup → index

###  Script to do it
Perfect! Here’s a **complete  script** that will take your **SAM files** for father, mother, and child, and do all the postprocessing steps in the correct order:

* Convert SAM → BAM
* Filter for properly paired reads
* Queryname sort → fixmate
* Coordinate sort → mark duplicates
* Index final BAMs

You can save this as, for example, `process_trio.sh` and run it.

```
#!/bin/
# Process mapped reads for a trio: father, mother, child
# Usage:  process_trio.sh

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
```

---

### **How to use:**

1. Save the script:

```
nano process_trio.sh
# Paste the script and save (Ctrl+O, Enter, Ctrl+X)
```

2. Make it executable:

```
chmod +x process_trio.sh
```

3. Run it:

```
 process_trio.sh
```

---

✅ After it finishes, you’ll have **deduplicated, properly filtered, and indexed BAMs** in:

```
~/OneDrive/Om/exomesequencing/mapped/processed/
```

Files:

* `father.dedup.bam` (+ index)
* `mother.dedup.bam` (+ index)
* `child.dedup.bam` (+ index)

These are ready for **variant calling** with FreeBayes or any other tool.

---

If you want, I can also **modify this script to clean up intermediate files** (like filtered and fixmate BAMs) so only the final deduplicated BAMs remain, saving disk space.

Do you want me to do that?

# VCF 
Perfect! Now that your **BAM files are fully processed**, we can move on to **variant calling with FreeBayes**. I’ll explain how to do this on your Linux system without Galaxy, using your deduplicated BAMs.

We have:

```
father.dedup.bam
mother.dedup.bam
child.dedup.bam
```

and your reference genome:

```
hg19_chr8.fa
```

---

## **Step 1: Install FreeBayes (if not already installed)**

```
sudo apt update
sudo apt install freebayes
```

Or, for the latest version:

```
conda install -c bioconda freebayes
```

---

## **Step 2: Run FreeBayes for multisample variant calling**

FreeBayes can accept multiple BAM files and output a single **VCF** file with genotypes for all samples.

```
freebayes -f ~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
father.dedup.bam mother.dedup.bam child.dedup.bam \
> trio_variants.vcf
```

**Explanation of flags:**

* `-f <reference.fa>` → the reference genome
* List of BAMs → the samples you want to call variants for (FreeBayes automatically merges them into a multisample VCF)
* Output is redirected to `trio_variants.vcf`

---

### **Optional parameters for simple diploid calling** (matching Galaxy “Simple diploid calling”):

```
freebayes -f hg19_chr8.fa \
--ploidy 2 \
--min-alternate-count 1 \
--min-alternate-fraction 0.2 \
father.dedup.bam mother.dedup.bam child.dedup.bam \
> trio_variants.vcf
```

* `--ploidy 2` → human diploid
* `--min-alternate-count 1` → call variants with at least 1 supporting read
* `--min-alternate-fraction 0.2` → at least 20% of reads must support the alternate allele

These are default-ish parameters for simple variant calling.

---

## **Step 3: Inspect the VCF**

The VCF is a **tab-delimited text file** that contains:

* Chromosome (`CHROM`)
* Position (`POS`)
* Reference allele (`REF`)
* Alternate allele(s) (`ALT`)
* Quality and filter info
* Genotypes for each sample (`GT` in the FORMAT column)

You can view it with:

```
less -S trio_variants.vcf
```

or

```
head trio_variants.vcf
```

For prettier inspection:

```
bcftools view trio_variants.vcf | less -S
```

---

### ✅ **Outcome**

* You now have a **multisample VCF** with all variants detected in father, mother, and child.
* Each sample’s genotype at each variant position is included.

---

If you want, I can provide a **full Linux script** that:

1. Takes your three deduplicated BAMs
2. Runs FreeBayes with diploid calling
3. Outputs a single VCF
4. Optionally compresses (`bgzip`) and indexes (`tabix`) the VCF for easier downstream analysis

This will mimic what Galaxy does automatically.

# Post Processing
Perfect! Now we’re moving to **post-processing your FreeBayes VCF** so it’s compatible with downstream annotation tools like SnpEff, GEMINI, or bcftools. On Linux, the equivalent of Galaxy’s `bcftools norm` workflow can be done directly with `bcftools`.

You already have:

```
trio_variants.vcf
hg19_chr8.fa
```

---

## **Step 1: Split multiallelic records & left-align indels**

We use `bcftools norm` with the following options:

```
bcftools norm -f ~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
-m -both \
--check-ref w \
-o trio_variants.norm.vcf \
trio_variants.vcf
```

### **Explanation of flags:**

| Flag            | Meaning                                                                             |
| --------------- | ----------------------------------------------------------------------------------- |
| `-f <ref.fa>`   | Reference genome used for normalization                                             |
| `-m -both`      | Split multiallelic sites into separate biallelic records (for both SNPs and indels) |
| `--check-ref w` | Ignore cases where REF does not match reference (`-w` in Galaxy)                    |
| `-o <output>`   | Output file                                                                         |
| Input VCF       | Your FreeBayes VCF                                                                  |

> `bcftools norm` automatically **left-aligns indels** and normalizes them according to the reference.

---

## **Step 2: Optional — view summary of changes**

You can check how many variants were split, realigned, or skipped using:

```
bcftools view trio_variants.norm.vcf | head
```

Or for a **detailed count**:

```
bcftools stats trio_variants.norm.vcf > trio_variants.norm.stats.txt
less trio_variants.norm.stats.txt
```

* `bcftools stats` reports total variant counts, SNPs vs indels, etc.

---

## **Step 3: Optional — compress & index normalized VCF**

```
bgzip trio_variants.norm.vcf
tabix -p vcf trio_variants.norm.vcf.gz
```

* This creates `trio_variants.norm.vcf.gz` and an index `trio_variants.norm.vcf.gz.tbi`
* Required for fast querying or tools like GEMINI.

---

✅ **Outcome:**

* `trio_variants.norm.vcf` → normalized, left-aligned, biallelic-split VCF ready for annotation
* Optionally compressed/indexed for efficient storage: `trio_variants.norm.vcf.gz`

---

If you want, I can write a **one-shot Linux script** that:

1. Takes your `trio_variants.vcf`
2. Runs `bcftools norm` to split & normalize
3. Compresses & indexes the result

This mimics the Galaxy workflow exactly.

Do you want me to do that?
## VCF Dependency issue
Ah, that’s expected on **modern Ubuntu versions** — `libssl1.0.0` is obsolete and no longer available in the official repositories. That’s why the system can’t find it.

The **best solution** here is to use **Conda**, which will provide a self-contained `bcftools` that doesn’t rely on the old library. This avoids all library conflicts.

---

### **Steps to fix with Conda**

1. **Create a Conda environment with bcftools:**

```
conda create -n bioinfo bcftools=1.17 samtools=1.17 -c bioconda -c conda-forge
```

2. **Activate the environment:**

```
conda activate bioinfo
```

3. **Check bcftools is working:**

```
bcftools --version
```

You should see something like:

```
bcftools 1.17
```

4. **Run your normalization command inside the Conda environment:**

```
bcftools norm -f ~/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
-m -both \
--check-ref w \
-o trio_variants.norm.vcf \
trio_variants.vcf
```

---

✅ After this, you’ll have a **normalized, left-aligned, biallelic-split VCF**, ready for annotation and downstream analysis, without any OpenSSL errors.

If you want, I can now write a **complete, ready-to-run Linux script** that does **postprocessed BAM → FreeBayes → bcftools norm → compress & index**, all in one shot using Conda.

Do you want me to do that?
Perfect! That means your reference FASTA is already indexed. We can now fix the VCF header and run `bcftools norm`. Using your paths, here’s a **ready-to-run sequence**:

---

### **Step 1: Get chr8 length from the index**

```
cut -f1,2 /home/thirdsem/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa.fai
```

* Example output:

```
chr8    146364022
```

---

### **Step 2: Create a contig header file**

```
echo "##contig=<ID=chr8,length=146364022>" > /home/thirdsem/OneDrive/Om/exomesequencing/mapped/contig_header.txt
```

---

### **Step 3: Prepend contig line to your VCF**

```
bcftools reheader -h <(cat /home/thirdsem/OneDrive/Om/exomesequencing/mapped/contig_header.txt <(grep '^#' /home/thirdsem/OneDrive/Om/exomesequencing/mapped/trio_variants.vcf)) \
-o /home/thirdsem/OneDrive/Om/exomesequencing/mapped/trio_variants.headered.vcf \
/home/thirdsem/OneDrive/Om/exomesequencing/mapped/trio_variants.vcf
```

---

### **Step 4: Normalize the VCF**

```
bcftools norm -f /home/thirdsem/OneDrive/Om/exomesequencing/reference_genome/hg19_chr8.fa \
-m -both \
--check-ref w \
-o /home/thirdsem/OneDrive/Om/exomesequencing/mapped/trio_variants.norm.vcf \
/home/thirdsem/OneDrive/Om/exomesequencing/mapped/trio_variants.headered.vcf
```

---

✅ After this, `trio_variants.norm.vcf` will be **left-aligned, normalized, and ready for annotation**.

If you want, I can also give the **final one-liner script** that does everything from header fix → normalization → compression → tabix indexing in one go. It’s very handy for your workflow.


