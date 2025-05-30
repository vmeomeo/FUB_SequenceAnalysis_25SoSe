rule iqtree:
    input:
        alignment = f"{config['output_dir_path']}/alignment/{{sample}}_{{read_number}}_aligned.fasta"
    output:
        treefile = f"{config['output_dir_path']}/phylogeny/{{sample}}_{{read_number}}/tree.treefile"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/phylogeny/iqtree_{{sample}}_{{read_number}}.log"
        # lambda wc : f"{config['output_dir_path']}/phylogeny/iqtree_{wc.sample}_{wc.read_number}.log"
    conda:
        "../env/iqtree.yaml"
    params:
        # prefix = lambda wc : f"{config['output_dir_path']}/phylogeny/{wc.sample}_{wc.read_number}/tree"
        prefix = f"{config['output_dir_path']}/phylogeny/{{sample}}_{{read_number}}/tree"
    shell:
        """
        mkdir -p $(dirname {log})
        iqtree -s {input.alignment} -nt {threads} -pre {params.prefix} > {log} 2>&1
        """