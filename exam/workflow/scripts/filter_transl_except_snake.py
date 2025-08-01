from Bio import SeqIO
import re

gff_in = snakemake.input.gff
faa_in = snakemake.input.faa
ffn_in = snakemake.input.ffn

gff_out = snakemake.output.gff
faa_out = snakemake.output.faa
ffn_out = snakemake.output.ffn

def get_invalid_ids(gff_file):
    invalid_ids = set()
    with open(gff_file) as f:
        for line in f:
            if "transl_except" in line:
                match = re.search(r"ID=([^;]+)", line)
                if match:
                    invalid_ids.add(match.group(1))
    return invalid_ids

def filter_fasta(input_file, output_file, invalid_ids):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record_id = record.id.split()[0]
            if record_id not in invalid_ids:
                SeqIO.write(record, out_handle, "fasta")

def filter_gff(input_file, output_file, invalid_ids):
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            match = re.search(r"ID=([^;]+)", line)
            if match and match.group(1) in invalid_ids:
                continue
            fout.write(line)

invalid_ids = get_invalid_ids(gff_in)
print(f"{gff_in} → {len(invalid_ids)} invalid gene(s)")
filter_gff(gff_in, gff_out, invalid_ids)
filter_fasta(faa_in, faa_out, invalid_ids)
filter_fasta(ffn_in, ffn_out, invalid_ids)
