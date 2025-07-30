rule panaroo:
    input:
        gffs = expand(
            f"{config['output_dir_path']}/annotation/{{sample}}/assembly.gff3",
            sample=samples.index #get_included_samples()
        )
    output:
        directory(f"{config['output_dir_path']}/pangenome/panaroo")
    threads: workflow.cores
    conda:
        "../env/panaroo.yaml"
    shell:
        """
        panaroo -i {input.gffs} \
                -o {output} \
                -t {threads} \
                --clean-mode strict \
                --align
        """
