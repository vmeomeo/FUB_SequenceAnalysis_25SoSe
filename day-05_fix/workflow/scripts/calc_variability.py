#!/usr/bin/env python3

import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from pathlib import Path

# Arguments
alignment_file = sys.argv[1]
window_size = int(sys.argv[2])
output_txt = sys.argv[3]
output_png = sys.argv[4]
step_size = int(sys.argv[5]) if len(sys.argv) > 5 else 1

# Read alignment
aln = AlignIO.read(alignment_file, "fasta")
length = aln.get_alignment_length()

# π calculation per column
def pi(column):
    pairs = list(itertools.combinations(column, 2))
    diffs = sum(1 for a, b in pairs if a != b and a != '-' and b != '-')
    return diffs / len(pairs) if pairs else 0

# Calculate π per column
pi_per_column = [pi(aln[:, i]) for i in range(length)]

# Rolling average with step size
windows = []
pis = []

for start in range(0, length - window_size + 1, step_size):
    end = start + window_size
    window_pi = np.mean(pi_per_column[start:end])
    windows.append(f"{start+1}-{end}")
    pis.append(window_pi)

# Write raw data to txt
with open(output_txt, "w") as f:
    f.write("Window\tPi\n")
    for w, v in zip(windows, pis):
        f.write(f"{w}\t{v:.6f}\n")

# Define thresholds
pi_array = np.array(pis)
low_threshold = np.percentile(pi_array, 10)
high_threshold = np.percentile(pi_array, 90)

low_indices = [i for i, v in enumerate(pis) if v <= low_threshold]
high_indices = [i for i, v in enumerate(pis) if v >= high_threshold]

# Write summary file
summary_file = Path(output_txt).with_suffix(".summary.txt")
with open(summary_file, "w") as f:
    f.write(f"Low variability windows (π ≤ {low_threshold:.4f}):\n")
    for i in low_indices:
        f.write(f"{windows[i]}\tπ = {pis[i]:.6f}\n")
    f.write(f"\nHigh variability windows (π ≥ {high_threshold:.4f}):\n")
    for i in high_indices:
        f.write(f"{windows[i]}\tπ = {pis[i]:.6f}\n")

# Plotting
plt.figure(figsize=(12, 5))
x = np.arange(len(pis))
colors = ['red' if i in high_indices else 'green' if i in low_indices else 'blue' for i in range(len(pis))]
plt.scatter(x, pis, c=colors, s=40)
plt.plot(x, pis, color='lightgray', linewidth=1)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label=f'High variability (top 10%, π ≥ {high_threshold:.4f})',
           markerfacecolor='red', markersize=10),
    Line2D([0], [0], marker='o', color='w', label=f'Low variability (bottom 10%, π ≤ {low_threshold:.4f})',
           markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Mid-variability (middle 80%)',
           markerfacecolor='blue', markersize=10)
]
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), framealpha=0.7)
plt.subplots_adjust(right=0.8)

tick_step = max(1, len(windows) // 20)
xtick_positions = np.arange(0, len(windows), tick_step)
xtick_labels = [windows[i] for i in xtick_positions]
plt.xticks(xtick_positions, xtick_labels, rotation=90)

plt.xlabel("Genome window (bp)")
plt.ylabel("Average pairwise nucleotide diversity (π)")
plt.title("Smoothed sequence variability across genome")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(output_png)
