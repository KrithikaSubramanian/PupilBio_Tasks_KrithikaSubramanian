# **Task 1: Data Handling and Statistical Analysis**
This repository contains Python scripts for the statistical analysis of phased methylation patterns (PMPs). Task 1 is divided into three steps, with corresponding Python scripts provided.

---

**Step 1: Coverage Analysis**
This script calculates the number of samples and tissues, along with the median and coefficient of variation (CV) for total CpG coverage within each tissue type.

```bash
python CoverageAnalysis.py
```

**Step 2: Identify Biomarkers**
This script identifies significant phased methylation patterns (PMPs) using the Chi-Square Test. It outputs:

A CSV file (Significant_PMPs_Data.csv) containing significant PMPs.
A summary of the unique PMP count and their average reads.
A bar plot (PMP_distributions.png) illustrating the distribution of total reads for significant PMPs.

```bash
python BiomarkerIdentifications.py
```

**Step 3: Calculate VRF**
This script calculates the Variant Read Fraction (VRF) for PMPs and generates two CSV files:
All_PMP_Mean_VRF.csv: Contains mean VRF for all PMPs.
Significant_PMP_Mean_VRF.csv: Contains mean VRF for significant PMPs.

```bash
python Calculate_VRF.py
```

# **Task 2: NGS Data Analysis**
This repository contains scripts for parallel processing of NGS data, somatic mutation calling, and statistical summarization. The workflow is designed to analyze `.fastq.gz` files using tools like **FASTQC**, **BWA**, **SAMTOOLS**, and **bcftools**, with the output summarized in an Excel file.

---

## **Table of Contents**
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Workflow](#workflow)
4. [Post-Analysis](#post-analysis)
5. [Repository Structure](#repository-structure)
6. [How to Run](#how-to-run)
7. [Contributing](#contributing)
8. [License](#license)

---

## **Overview**

The workflow consists of five main steps:
1. **FastQC**: Perform quality control checks on raw sequencing data.
2. **BWA Indexing and Alignment**: Index the reference genome and map reads to the reference.
3. **BAM Processing**: Convert SAM files to BAM format, sort BAM files, and create indexes.
4. **Somatic Mutation Calling**: Identify somatic mutations using GATK or similar tools.
5. **VCF Processing**: Extract key mutation information from VCF files using `bcftools`.

---

## **Prerequisites**

Ensure the following tools and libraries are installed:

### **Software Requirements**
- **Python 3.x**
- **FASTQC**
- **BWA**
- **SAMTOOLS**
- **bcftools**
- **MPI** (for parallel processing)

### **Python Libraries**
Install the required Python libraries:
```bash
pip install pandas numpy openpyxl
```

## **Workflow**

Execution Steps
Steps 1–3: Run in parallel for both samples using the following command:

```bash
mpirun -np 2 python Task2_NGS_script_Krithika.py
```

Steps 4–5: Run sequentially without MPI:
```bash
python Task2_NGS_script_Krithika.py
```

**Steps Explained**
- **FastQC:** Perform quality control checks on raw .fastq.gz files to assess read quality.
- **BWA Index and BWA MEM:**Index the reference genome. Map sequencing reads to the reference genome to generate SAM files.
- **BAM Processing:** Convert SAM to BAM format.Sort BAM files and index them for downstream analysis.
- **Somatic Mutation Calling:** Identify somatic mutations using GATK-Mutect2 or similar tools.
- **VCF Processing:** Extract relevant mutation information from the VCF file using bcftools.

## **Post-Analysis**
To calculate summary statistics for the mutations:
Execute the statistics calculation script:

```bash
python CalculateStatistics.py
```

This will generate an Excel file named variant_analysis.xlsx, which contains a detailed statistical summary.

## **How to Run**

**Set up your files:**
Place your raw .fastq.gz files in the input/ directory.
Place the reference genome .fasta file in the reference/ directory.

**Run the workflow:**
For parallel processing of the first three steps:
```bash
mpirun -np 2 python Task2_NGS_script_Krithika.py
```
For steps 4 and 5, execute sequentially:
```bash
python Task2_NGS_script_Krithika.py
```

Generate statistics:
Run the statistics script to create an Excel summary file:
```bash
python CalculateStatistics.py
```

