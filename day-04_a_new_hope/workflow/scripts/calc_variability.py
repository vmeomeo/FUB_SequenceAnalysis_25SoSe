#!/usr/bin/env python3
import sys
import itertools
import matplotlib.pyplot as plt
from Bio import AlignIO

alignment_file = sys.argv[1]
window_size = int(sys.argv[2])
output_txt = sys.argv[3]
output_png = sys.argv[4]

aln = AlignIO.read(alignment_file, "fasta")
length = aln.get_alignment_length()

def pi(column):
    pairs = list(itertools.combinations(column, 2))
    diffs = sum(1 for a, b in pairs if a != b and a != '-' and b != '-')
    return diffs / len(pairs) if pairs else 0

windows = []
pis = []

for start in range(0, length, window_size):
    end = min(start + window_size, length)
    if end - start == 0:
        continue
    window_cols = [aln[:, i] for i in range(start, end)]
    window_pi = sum(pi(col) for col in window_cols) / (end - start)
    windows.append(f"{start+1}-{end}")
    pis.append(window_pi)

with open(output_txt, "w") as f:
    f.write("Window\tPi\n")
    for w, v in zip(windows, pis):
        f.write(f"{w}\t{v:.6f}\n")

plt.figure(figsize=(10,4))
plt.plot(range(len(pis)), pis, marker='o')
plt.xticks(range(len(windows)), windows, rotation=90)
plt.xlabel("Genome window (bp)")
plt.ylabel("Average pairwise nucleotide diversity (π)")
plt.title("Sequence variability across genome")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(output_png)
