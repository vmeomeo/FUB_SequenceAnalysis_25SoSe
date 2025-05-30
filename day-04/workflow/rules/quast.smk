rule quast:
    input: 
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}_{{read_number}}/contigs.fasta"
    output:
        assembly_qc = directory(f"{config['output_dir_path']}/assembly_qc/{{sample}}_{{read_number}}")
    conda:
        "../env/quast.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/assembly_qc/logs/quast_{{sample}}_{{read_number}}.log"
    shell:
        """
        quast.py {input.assembly} -o {output.assembly_qc} -t {threads} > {log} 2>&1
        """