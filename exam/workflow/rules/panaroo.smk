# rule filter_gffs:
#     input:
#         gff = f"{config['output_dir_path']}/annotation/{{sample}}/assembly.gff3",
#         faa = f"{config['output_dir_path']}/annotation/{{sample}}/assembly.faa",
#         ffn = f"{config['output_dir_path']}/annotation/{{sample}}/assembly.ffn"
#     output:
#         gff = temp("results/annotation/{sample}/assembly.filtered.gff3"),
#         faa = temp("results/annotation/{sample}/assembly.filtered.faa"),
#         ffn = temp("results/annotation/{sample}/assembly.filtered.ffn"),
#     conda:
#         "../env/bioseq.yaml"
#     script:
#         "workflow/scripts/filter_transl_except_snake.py"

rule panaroo:
    input:
        # gffs = expand("results/annotation/{sample}/assembly.filtered.gff3", sample=samples.index)
        # gffs = expand(f"{config['output_dir_path']}/annotation/{{sample}}/assembly.gff3", sample=samples.index)
        # Test
        gffs = expand(f"{config['output_dir_path']}/panaroo_input/{{sample}}.gff3", sample=samples.index)
    output:
        # directory(f"{config['output_dir_path']}/pangenome/panaroo")
        # Test
        directory(f"{config['output_dir_path']}/pangenome/panaroo_test")
    # log:
        # f"{config['output_dir_path']}/logs/panaroo/{{sample}}.log"
    conda: "../env/panaroo.yaml"
    threads: workflow.cores
    shell:
        """
        panaroo -i {input.gffs} \
                -o {output} \
                -t {threads} \
                --clean-mode strict \
                --merge_paralogs \
                --aligner mafft \
                --alignment pan \
                --remove-invalid-genes
        """
        #               --len_dif_percent 20 \ The culprit of the error was this line
        #  > {log} 2>&1
