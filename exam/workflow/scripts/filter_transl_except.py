from Bio import SeqIO
import os
import re
import sys

samples = ["sample1", "sample2", "sample3", "sample4", "sample5"]
annotation_path = "results/annotation"

def get_invalid_ids(gff_file):
    """Return a set of gene IDs that have transl_except annotations."""
    invalid_ids = set()
    with open(gff_file) as f:
        for line in f:
            if "transl_except" in line:
                match = re.search(r"ID=([^;]+)", line)
                if match:
                    invalid_ids.add(match.group(1))
    return invalid_ids

def filter_fasta(input_file, output_file, invalid_ids):
    """Write only sequences whose IDs are not in invalid_ids."""
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record_id = record.id.split()[0]
            if record_id not in invalid_ids:
                SeqIO.write(record, out_handle, "fasta")

def filter_gff(input_file, output_file, invalid_ids):
    """Write lines from the GFF file except those matching IDs in invalid_ids."""
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            match = re.search(r"ID=([^;]+)", line)
            if match and match.group(1) in invalid_ids:
                continue
            fout.write(line)

for sample in samples:
    gff = os.path.join(annotation_path, sample, "assembly.gff3")
    faa = os.path.join(annotation_path, sample, "assembly.faa")
    ffn = os.path.join(annotation_path, sample, "assembly.ffn")

    filtered_gff = gff.replace(".gff3", ".filtered.gff3")
    filtered_faa = faa.replace(".faa", ".filtered.faa")
    filtered_ffn = ffn.replace(".ffn", ".filtered.ffn")

    bad_ids = get_invalid_ids(gff)
    print(f"{sample}: {len(bad_ids)} invalid gene(s)")

    filter_gff(gff, filtered_gff, bad_ids)
    filter_fasta(faa, filtered_faa, bad_ids)
    filter_fasta(ffn, filtered_ffn, bad_ids)
