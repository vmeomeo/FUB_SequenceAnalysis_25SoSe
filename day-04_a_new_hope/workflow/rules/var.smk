rule var:
    input:
        alignment = f"{config['output_dir_path']}/alignment/all_samples_aligned.fasta"
    output:
        variability_txt = f"{config['output_dir_path']}/variability/all_samples_windows_var.txt",
        plot_png = f"{config['output_dir_path']}/variability/all_samples_var_plot.png"
    params:
        window_size = config['window_size']
    conda:
        "../env/var.yaml"
    log:
        f"{config['output_dir_path']}/variability/all_samples_variability.log"
    shell:
        """
        mkdir -p $(dirname {output.variability_txt})

        python -c "
import sys
import itertools
import matplotlib.pyplot as plt
from Bio import AlignIO

window_size = {params.window_size}
aln = AlignIO.read('{input.alignment}', 'fasta')
length = aln.get_alignment_length()

def pi(column):
    pairs = list(itertools.combinations(column, 2))
    diffs = sum(1 for a, b in pairs if a != b and a != '-' and b != '-')
    return diffs / len(pairs) if pairs else 0

windows = []
pis = []

for start in range(0, length, window_size):
    end = min(start + window_size, length)
    window_cols = [aln[:, i] for i in range(start, end)]
    window_pi = sum(pi(col) for col in window_cols) / (end - start)
    windows.append(f'{{start+1}}-{{end}}')
    pis.append(window_pi)

with open('{output.variability_txt}', 'w') as f:
    f.write('Window\\tPi\\n')
    for w, v in zip(windows, pis):
        f.write(f'{{w}}\\t{{v:.6f}}\\n')

plt.figure(figsize=(10,4))
plt.plot(range(len(pis)), pis, marker='o')
plt.xticks(range(len(windows)), windows, rotation=90)
plt.xlabel('Genome window (bp)')
plt.ylabel('Average pairwise nucleotide diversity (π)')
plt.title('Sequence variability across genome')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('{output.plot_png}')
" > {log} 2>&1
        """
