from Bio import SeqIO
from pathlib import Path

# Liste des gènes posant problème
invalid_genes = {
    "sample1": "ILKLBK_01339",
    "sample2": "IHEPIH_03111",
    "sample3": "NKHLMP_02384",
    "sample4": "DHNBJF_02282",
    "sample5": "BPEKIH_02216"
}

base_path = Path("results/annotation")

for sample, gene_id in invalid_genes.items():
    ffn_file = base_path / sample / "assembly.ffn"
    print(f"Checking {sample} - {gene_id}")
    found = False
    for record in SeqIO.parse(ffn_file, "fasta"):
        if gene_id in record.id:
            seq = str(record.seq)
            found = True
            print(f"Length: {len(seq)}")
            print(f"Start codon: {seq[:3]}")
            print(f"Stop codon: {seq[-3:]}")
            print(f"Divisible by 3: {len(seq) % 3 == 0}")
            break
    if not found:
        print(f"{gene_id} not found in {ffn_file}")
    print("-" * 40)
