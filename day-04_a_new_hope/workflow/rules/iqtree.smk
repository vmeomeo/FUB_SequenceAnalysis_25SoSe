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
        prefix = f"{config['output_dir_path']}/phylogeny/all_samples",
        threads_max = 32
    shell:
        """
        mkdir -p $(dirname {log})
        iqtree -s {input.alignment} -nt {threads} -ntmax {params.threads_max} -pre {params.prefix} > {log} 2>&1
        """
