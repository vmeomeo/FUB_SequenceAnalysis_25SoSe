rule spades:
    input:
        r1 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_1.fastq",
        r2 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_2.fastq"
    output:
        assembly_dir = directory(f"{config['output_dir_path']}/assembly/{{sample}}"),
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}/contigs.fasta"
    conda:
        "../env/spades.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/assembly/logs/spades_{{sample}}.log"
    shell:
        """
        spades.py -1 {input.r1} -2 {input.r2} -o {output.assembly_dir} --isolate -t {threads} > {log} 2>&1
        """
