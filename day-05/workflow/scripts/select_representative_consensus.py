#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO

input_files = sys.argv[1:-1]
output_file = sys.argv[-1]

os.makedirs(os.path.dirname(output_file), exist_ok=True)

representative_seqs = {}

for filepath in input_files:
    filename = os.path.basename(filepath)
    # Remove the suffix to isolate sample and read number
    basename = filename.replace(".consensus.fasta", "")  # e.g. sample_22_1
    parts = basename.split("_")
    sample = "_".join(parts[:-1])   # e.g. "sample_22"
    read_number = parts[-1]         # e.g. "1"
    full_sample_id = f"{sample}_{read_number}"

    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        continue

    # Select longest sequence from the fasta file
    longest = max(records, key=lambda r: len(r.seq))
    longest.id = full_sample_id
    longest.description = ""

    # Keep only the longest representative sequence per sample
    if sample not in representative_seqs or len(longest.seq) > len(representative_seqs[sample].seq):
        representative_seqs[sample] = longest

with open(output_file, "w") as out_f:
    SeqIO.write(representative_seqs.values(), out_f, "fasta")
