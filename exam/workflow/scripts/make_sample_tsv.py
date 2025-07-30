import os
import pandas as pd
from pathlib import Path

# Set this to your directory containing FASTQ files (relative to Snakemake file)
fastq_dir_set1 = Path("data/SET1")
fastq_dir_set2 = Path("data/SET2")

# Find all _1.fastq.gz and _2.fastq.gz files
fq1_files_set2 = sorted([f for f in fastq_dir_set2.glob("*_1.fastq.gz")])
fq2_files_set2 = sorted([f for f in fastq_dir_set2.glob("*_2.fastq.gz")])

# Match by sample name prefix
samples = []
for fq1 in fq1_files_set2:
    sample_name = fq1.name.replace("_1.fastq.gz", "")
    fq2 = fastq_dir_set2 / f"{sample_name}_2.fastq.gz"
    if fq2.exists():
        samples.append({
            "sample": sample_name,
            "fq1": str(fq1),
            "fq2": str(fq2)
        })

# Create and save samples.tsv
df = pd.DataFrame(samples)
df.to_csv("config/samples_set2.tsv", sep="\t", index=False)
print("samples_set2.tsv created with", len(df), "samples.")

# for set1
# Find all _1.fastq.gz and _2.fastq.gz files
fq1_files_set1 = sorted([f for f in fastq_dir_set1.glob("*_1.fastq.gz")])
fq2_files_set1 = sorted([f for f in fastq_dir_set1.glob("*_2.fastq.gz")])
fq3_files_set1 = sorted([f for f in fastq_dir_set1.glob("*_sequencing.fastq.gz")])

# Match by sample name prefix
samples = []
for fq1 in fq1_files_set1:
    sample_name = fq1.name.replace("_Illumina_MiSeq_paired_end_sequencing_1.fastq.gz", "")
    fq2 = fastq_dir_set1 / f"{sample_name}_Illumina_MiSeq_paired_end_sequencing_2.fastq.gz"
    fq3 = fastq_dir_set1 / f"{sample_name}_MinION_sequencing.fastq.gz"
    if (fq2.exists() and fq3.exists()):
        samples.append({
            "sample": sample_name,
            "fq1": str(fq1),
            "fq2": str(fq2),
            "fq3": str(fq3)
        })

# Create and save samples.tsv
df = pd.DataFrame(samples)
df.to_csv("config/samples_set1.tsv", sep="\t", index=False)
print("samples_set1.tsv created with", len(df), "samples.")