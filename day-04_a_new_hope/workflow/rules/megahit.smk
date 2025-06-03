rule megahit:
    input:
        r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_1.fastq",
        r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_2.fastq"
    output:
        assembly_dir = directory(f"{config['output_dir_path']}/assembly/{{sample}}")
    conda:
        "../env/megahit.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/assembly/logs/megahit_{{sample}}.log"
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -o {output.assembly_dir} -t {threads} > {log} 2>&1
        """
