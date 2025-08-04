# TL;DR
### How to run snakemake?

1. Go to directory outside 'workflow' 1 level

2. Dry-run snakemake to check errors
```snakemake --use-conda --cores 16 -n```

3. If there are no errors, run (you can change the number of cores)
```snakemake --use-conda --cores 16```

This Snakemake pipeline performs quality control, hybrid or short-read-only assembly, annotation, and genomic feature prediction for bacterial samples (e.g. *Klebsiella pneumoniae*). It supports both short reads (Illumina) and long reads (Nanopore) depending on your dataset.

---
# Detail on the workflow
## 📁 Directory Structure

```
project/
├── config/
│   ├── samples_set1.tsv       # Created by script (for long+short reads)
│   └── samples_set2.tsv       # Created by script (for short reads only)
├── data/
│   ├── SET1/                  # Contains raw FASTQ files for Set 1
│   └── SET2/                  # Contains raw FASTQ files for Set 2
├── resources/                 # Databases (e.g. Bakta, BUSCO, etc.)
├── workflow/  
│   ├── env                    # Conda environment configurations for each rule
│   ├── rules                  # Snakemake rules
│   ├── scripts/               # Helper scripts (see below)
│   └── Snakefile              # Main snakemake file
├── results/                   # Output directory (auto-created)
└── README.md                  # Tutorial on how to use this workflow
```

---

## 🛠️ Configuration

Edit `config.yaml` to adjust:

- Input mode: `samples: "config/samples_set1.tsv"` or `samples_set2.tsv`
- Read type: `short_reads` and `long_reads` (boolean)
- Enable/disable tools: `run_mlst`, `run_resistance`, etc.
- Paths to reference databases: `bakta_db_path`, `busco_db_path`, ...

---

## 🧾 Creating Sample Sheets

Use the script [`scripts/make_sample_tsv.py`](workflow/scripts/make_sample_tsv.py) to automatically generate TSV files for your datasets.

Run:

```bash
python scripts/make_sample_tsv.py
```

This generates:

- `config/samples_set1.tsv` for samples with 3 files (fq1, fq2, fq3)
- `config/samples_set2.tsv` for short-read-only samples

Use can customize for your own dataset, by placing your data folder in the `data` directory, and edit the path and variables using the above python script as template.

---

## 🚀 How to Run the Pipeline

1. **Navigate to project root directory** (one level above `workflow/`):

```bash
cd path/to/project/
```

2. **Dry-run to check for errors:**

```bash
snakemake --use-conda --cores 16 -n
```

3. **Execute workflow:**

```bash
snakemake --use-conda --cores 16
```

(You can change `--cores` as needed.)

---

## 📚 Notes

- **Blacklist**: Optional. If you want to exclude samples, create `config/sample_blacklist.tsv` with a `sample` column.
- If `long_reads: false`, the pipeline will skip long-read steps (e.g., Unicycler `-l`, NanoPlot QC).
- BUSCO/Bakta DBs will be auto-downloaded if `*_db_download: true`.