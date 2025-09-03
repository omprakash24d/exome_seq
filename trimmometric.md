
# 📄 **Documentation: Installing and Running Trimmomatic for Exome Sequencing (Paired-End Data)**

---

## 📌 Overview

**Goal:**
To perform quality trimming and adapter removal on paired-end FASTQ files using **Trimmomatic**, as part of an **exome sequencing pipeline** for three samples: `father`, `mother`, and `proband`.

---

## 1️⃣ Installation of Trimmomatic via Conda

### ✅ What?

Trimmomatic is a flexible read trimming tool for Illumina NGS data. It removes adapters, trims low-quality ends, and filters short reads.

### 💡 Why Conda?

Using Conda (or Mamba) ensures clean installation and manages dependencies easily — especially useful in bioinformatics workflows.

### 💻 How?

We installed it using Bioconda in a dedicated Conda environment called `exome-seq`.

```bash
conda create -n exome-seq
conda activate exome-seq
conda install -c bioconda trimmomatic
```

### 📍 Where is it installed?

Trimmomatic and its adapter files are located inside the Conda environment:

```
/home/thirdsem/miniconda3/envs/exome-seq/share/trimmomatic-0.39-2/
```

---

## 2️⃣ Input FASTQ Files

Stored in:

```bash
~/OneDrive/Om/exomesequencing/raw_data
```

Files (paired-end format):

```bash
father_R1.fq.gz   father_R2.fq.gz
mother_R1.fq.gz   mother_R2.fq.gz
proband_R1.fq.gz  proband_R2.fq.gz
```

---

## 3️⃣ Output Directory

Created to store the cleaned FASTQ files:

```bash
~/OneDrive/Om/exomesequencing/trimmed_data/
```

---

## 4️⃣ Trimming Command Breakdown

### ✂️ Tool Used:

```bash
trimmomatic PE
```

### 🔍 What each parameter does:

| Parameter            | Purpose                                            |
| -------------------- | -------------------------------------------------- |
| `PE`                 | Paired-end mode                                    |
| `-threads 4`         | Use 4 CPU cores for faster processing              |
| `-phred33`           | Force interpretation of quality scores as Phred+33 |
| `ILLUMINACLIP`       | Remove Illumina adapters using adapter file        |
| `LEADING:3`          | Remove leading bases with quality < 3              |
| `TRAILING:3`         | Remove trailing bases with quality < 3             |
| `SLIDINGWINDOW:4:15` | Cut when 4-base window average quality < 15        |
| `MINLEN:36`          | Drop reads shorter than 36 bases after trimming    |

### 🔬 Adapter file used:

```bash
/home/thirdsem/miniconda3/envs/exome-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
```

---

## 5️⃣ Full Trimming Commands (Used)

### 🧬 Father:

```bash
trimmomatic PE -threads 4 \
  father_R1.fq.gz father_R2.fq.gz \
  ../trimmed_data/father_R1_paired.fq.gz ../trimmed_data/father_R1_unpaired.fq.gz \
  ../trimmed_data/father_R2_paired.fq.gz ../trimmed_data/father_R2_unpaired.fq.gz \
  ILLUMINACLIP:/home/thirdsem/miniconda3/envs/exome-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 🧬 Mother:

```bash
trimmomatic PE -threads 4 \
  mother_R1.fq.gz mother_R2.fq.gz \
  ../trimmed_data/mother_R1_paired.fq.gz ../trimmed_data/mother_R1_unpaired.fq.gz \
  ../trimmed_data/mother_R2_paired.fq.gz ../trimmed_data/mother_R2_unpaired.fq.gz \
  ILLUMINACLIP:/home/thirdsem/miniconda3/envs/exome-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 🧬 Proband (Special Case):

The first run failed due to this error:

```
Error: Unable to detect quality encoding
```

#### ✅ Diagnosis:

* Valid FASTQ format ✅
* But many reads had **full `#` lines** (low-quality placeholders)
* Trimmomatic failed to auto-detect Phred encoding

#### ✅ Fix: Forced Phred+33

```bash
trimmomatic PE -threads 4 -phred33 \
  proband_R1.fq.gz proband_R2.fq.gz \
  ../trimmed_data/proband_R1_paired.fq.gz ../trimmed_data/proband_R1_unpaired.fq.gz \
  ../trimmed_data/proband_R2_paired.fq.gz ../trimmed_data/proband_R2_unpaired.fq.gz \
  ILLUMINACLIP:/home/thirdsem/miniconda3/envs/exome-seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

---

## ✅ Output Files (Per Sample)

Each sample produces:

* `<sample>_R1_paired.fq.gz` → High-quality, forward reads
* `<sample>_R1_unpaired.fq.gz` → Forward reads without a pair
* `<sample>_R2_paired.fq.gz` → High-quality, reverse reads
* `<sample>_R2_unpaired.fq.gz` → Reverse reads without a pair

All saved in:

```bash
~/OneDrive/Om/exomesequencing/trimmed_data/
```

---

## 🧪 Post-Trimming Check (Recommended)

After trimming, assess the quality using **FastQC**:

```bash
fastqc ../trimmed_data/*_paired.fq.gz
```

---

## 🧠 Tips for Reproducibility

* Use a bash script to automate trimming of multiple samples
* Keep a copy of your commands in a log file
* Store versions of tools (Trimmomatic v0.39, Conda env)
* Keep original and trimmed data clearly separated

---

## 📚 References

* [Trimmomatic GitHub](https://github.com/usadellab/Trimmomatic)
* [Bioconda Trimmomatic Recipe](https://bioconda.github.io/recipes/trimmomatic/README.html)
* FastQ format: [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format)

---

