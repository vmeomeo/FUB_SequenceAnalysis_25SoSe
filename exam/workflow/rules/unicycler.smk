rule unicycler:
    input:
        r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_r1.fastq.gz",
        r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_r2.fastq.gz"
    output:
        f"{config['output_dir_path']}/assembly/{{sample}}/assembly.fasta"
    log:
        f"{config['output_dir_path']}/logs/unicycler/{{sample}}.log"
    threads: workflow.cores
    conda:
        "../env/unicycler.yaml"
    params:
        outdir = config['output_dir_path']
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -o {params.outdir}/assembly/{wildcards.sample} "
        "-t {threads} &> {log}"