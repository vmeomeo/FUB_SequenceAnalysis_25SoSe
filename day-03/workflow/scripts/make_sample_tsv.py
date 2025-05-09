import os
import pandas as pd
from pathlib import Path

# Set this to your directory containing FASTQ files (relative to Snakemake file)
fastq_dir = Path("../data/fastq/small/subsampled_0.1")

# Find all _1.fastq.gz and _2.fastq.gz files
fq1_files = sorted([f for f in fastq_dir.glob("*_1.fastq.gz")])
fq2_files = sorted([f for f in fastq_dir.glob("*_2.fastq.gz")])

# Match by sample name prefix
samples = []
for fq1 in fq1_files:
    sample_name = fq1.name.replace("_1.fastq.gz", "")
    fq2 = fastq_dir / f"{sample_name}_2.fastq.gz"
    if fq2.exists():
        samples.append({
            "sample": sample_name,
            "fq1": str(fq1),
            "fq2": str(fq2)
        })

# Create and save samples.tsv
df = pd.DataFrame(samples)
df.to_csv("../config/samples.tsv", sep="\t", index=False)
print("samples.tsv created with", len(df), "samples.")
