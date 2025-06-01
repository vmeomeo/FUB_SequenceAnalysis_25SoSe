rule spades:
    input:
        cleaned = f"{config['output_dir_path']}/cleaned/{{sample}}_{{read_number}}.fastq"
    output:
        assembly_dir = directory(f"{config['output_dir_path']}/assembly/{{sample}}_{{read_number}}")
    conda:
        "../env/spades.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/assembly/logs/spades_{{sample}}_{{read_number}}.log"
    shell:
        """
        spades.py -s {input.cleaned} -o {output.assembly_dir} --isolate -t {threads} > {log} 2>&1
        """