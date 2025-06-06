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
        mkdir -p $(dirname {output.image})
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
