rule panaroo:
    input:
        gffs = expand(f"{config['output_dir_path']}/annotation/{{sample}}/{{sample}}.gff3", sample=get_included_samples())
    output:
        dir = directory(f"{config['output_dir_path']}/pangenome/panaroo"),
        core = f"{config['output_dir_path']}/pangenome/panaroo/core_gene_alignment.aln"
    log:
        f"{config['output_dir_path']}/logs/panaroo/panaroo.log"
    conda: "../env/panaroo.yaml"
    threads: workflow.cores
    shell:
        """
        panaroo -i {input.gffs} \
                -o {output.dir} \
                -t {threads} \
                --clean-mode strict \
                --merge_paralogs \
                --aligner mafft \
                --alignment pan \
                --remove-invalid-genes \
                > {log} 2>&1
        """
