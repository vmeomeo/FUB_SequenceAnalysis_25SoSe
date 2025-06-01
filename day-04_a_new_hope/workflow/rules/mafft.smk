rule mafft:
    input:
        fasta = f"{config['output_dir_path']}/all_samples/all_samples.fasta"
    output:
        alignment = f"{config['output_dir_path']}/alignment/all_samples_aligned.fasta"
    conda:
        "../env/mafft.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/alignment/logs/mafft_all_samples.log"
    shell:
        """
        mkdir -p $(dirname {log})
        mafft --auto --thread {threads} {input.fasta} > {output.alignment} 2> {log}
        """
