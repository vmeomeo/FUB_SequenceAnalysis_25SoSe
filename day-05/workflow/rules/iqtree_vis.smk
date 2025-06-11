rule iqtree:
    input:
        alignment = f"{config['output_dir_path']}/alignment/all_samples_aligned.fasta"
    output:
        treefile = f"{config['output_dir_path']}/phylogeny/all_samples.treefile"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/phylogeny/iqtree_all_samples.log"
    conda:
        "../env/iqtree.yaml"
    params:
        prefix = f"{config['output_dir_path']}/phylogeny/all_samples"
    shell:
        """
        iqtree -s {input.alignment} -nt {threads} -pre {params.prefix} > {log} 2>&1
        """
rule visualize_tree:
    input:
        treefile = f"{config['output_dir_path']}/phylogeny/all_samples.treefile"
    output:
        image = f"{config['output_dir_path']}/phylogeny/all_samples.tree.png"
    conda:
        "../env/vis.yaml"
    log:
        f"{config['output_dir_path']}/phylogeny/ete3_all_samples.log"
    shell:
        """
        xvfb-run -a python -c "
import sys
from ete3 import Tree

t = Tree('{input.treefile}', format=0)

# Rename leaf names by removing '_1'
for leaf in t.iter_leaves():
    if leaf.name.endswith('_1'):
        leaf.name = leaf.name[:-2]

t.render('{output.image}', w=800, units='px')
" > {log} 2>&1
        """
