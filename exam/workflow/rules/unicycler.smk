rule unicycler:
    input:
        r1 = lambda wc: (
            f"{config['output_dir_path']}/decontaminated/{wc.sample}_r1.fastq"
            if config.get("contaminant_fasta", "")
            else f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz"
        ),
        r2 = lambda wc: (
            f"{config['output_dir_path']}/decontaminated/{wc.sample}_r2.fastq"
            if config.get("contaminant_fasta", "")
            else f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz"
        ),
        l = lambda wc: (
            f"{config['output_dir_path']}/filtered/{wc.sample}_filtered.fastq"
            if config.get("long_reads", False) else None
        )
    output:
        f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta"
    log:
        f"{config['output_dir_path']}/logs/unicycler/{{sample}}.log"
    threads: workflow.cores
    conda:
        "../env/unicycler.yaml"
    params:
        outdir = config['output_dir_path']
    shell:
        """
        mkdir -p {params.outdir}/assembly/{wildcards.sample}/temp
        unicycler -1 {input.r1} -2 {input.r2} \
            {"-l " + input.l if input.l else ""} \
            -o {params.outdir}/assembly/{wildcards.sample}/temp -t {threads} \
            &> {log}
        mv {params.outdir}/assembly/{wildcards.sample}/temp/assembly.fasta {output}
        """