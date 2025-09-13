rule iqtree:
    input:
        alignment = f"{config['output_dir_path']}/pangenome/panaroo/core_gene_alignment.aln"
    output:
        treefile = f"{config['output_dir_path']}/phylogeny/panaroo_tree.treefile"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/phylogeny/iqtree_panaroo.log"
    conda:
        "../env/iqtree.yaml"
    params:
        prefix = f"{config['output_dir_path']}/phylogeny/panaroo_tree"
    shell:
        """
        # Check if alignment has at least 4 sequences
        NUM_SEQS=$(grep -c '^>' {input.alignment} || wc -l < {input.alignment})
        if [ "$NUM_SEQS" -lt 4 ]; then
            echo "WARNING: Only $NUM_SEQS sequences detected. Skipping bootstrapping." > {log}
            iqtree2 -s {input.alignment} \
                    -nt {threads} \
                    -pre {params.prefix} \
                    -m MFP \
                    -n 0  # Disable bootstrapping for small datasets >> {log} 2>&1
        else
            iqtree2 -s {input.alignment} \
                    -nt {threads} \
                    -pre {params.prefix} \
                    -m MFP \
                    -B 1000 >> {log} 2>&1
        fi
        """

rule visualize_tree:
    input:
        treefile = f"{config['output_dir_path']}/phylogeny/panaroo_tree.treefile"
    output:
        image = f"{config['output_dir_path']}/phylogeny/panaroo_tree.png"
    conda:
        "../env/vis.yaml"
    log:
        f"{config['output_dir_path']}/phylogeny/ete3_panaroo.log"
    shell:
        """
        xvfb-run -a python -c "
from ete3 import Tree
t = Tree('{input.treefile}', format=0)
t.render('{output.image}', w=800, units='px')
" > {log} 2>&1
        """
