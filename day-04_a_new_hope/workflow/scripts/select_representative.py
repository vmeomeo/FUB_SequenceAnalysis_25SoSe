#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from collections import defaultdict

scaffold_files = sys.argv[1:-1]
output_file = sys.argv[-1]

os.makedirs(os.path.dirname(output_file), exist_ok=True)
representative_seqs = {}

for filepath in scaffold_files:
    filename = os.path.basename(filepath)
    parts = filename.split("_")
    sample = "_".join(parts[:-2])
    read_number = parts[-2]
    full_sample_id = f"{sample}_{read_number}"

    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        continue

    longest = max(records, key=lambda r: len(r.seq))
    longest.id = full_sample_id
    longest.description = ""

    if sample not in representative_seqs or len(longest.seq) > len(representative_seqs[sample].seq):
        representative_seqs[sample] = longest

with open(output_file, "w") as out_f:
    SeqIO.write(representative_seqs.values(), out_f, "fasta")
