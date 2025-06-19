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
    # Extract sample name from file name (e.g. "sample_22" from "sample_22.consensus.fasta")
    sample = filename.replace(".consensus.fasta", "")

    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        continue

    # Choose longest sequence from the file
    longest = max(records, key=lambda r: len(r.seq))
    longest.id = sample
    longest.description = ""

    representative_seqs[sample] = longest

# Write combined output
with open(output_file, "w") as out_f:
    SeqIO.write(representative_seqs.values(), out_f, "fasta")
