rule mafft:
    input: 
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}_{{read_number}}/contigs.fasta"
    output: 
        alignment = f"{config['output_dir_path']}/alignment/{{sample}}_{{read_number}}_aligned.fasta"
    conda:
        "../env/mafft.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/alignment/logs/mafft_{{sample}}_{{read_number}}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        mafft --auto --thread {threads} {input.assembly} > {output.alignment} 2> {log}
        """