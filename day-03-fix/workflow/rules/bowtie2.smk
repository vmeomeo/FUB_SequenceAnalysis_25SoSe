rule create_index:
    input: 
        lambda wc: f"{config['ref']}"
    output:
        expand(
            "{output_dir}/ref_idx/reference_index.{ext}",
            output_dir=config["output_dir_path"],
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    threads: workflow.cores
    conda:
        "../env/bowtie2.yaml"
    log:
        f"{config['output_dir_path']}/logs/create_index.log"
    shell:
        "bowtie2-build -f --threads {threads} {input} {config[output_dir_path]}/ref_idx/reference_index > {log} 2>&1"

rule align:
    input: 
        fq1 = f"{config['output_dir_path']}/trimmed/{{sample}}_R1_paired.fastq.gz",
        fq2 = f"{config['output_dir_path']}/trimmed/{{sample}}_R2_paired.fastq.gz",
        index = f"{config['output_dir_path']}/ref_idx/reference_index.1.bt2"  # dummy dependency
    output:
        sam = f"{config['output_dir_path']}/sam/{{sample}}.sam"
    threads: workflow.cores
    conda:
        "../env/bowtie2.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_bowtie2_align.log"
    params:
        index = f"{config['output_dir_path']}/ref_idx/reference_index",
        N = config["bowtie2"]["N"],
        L = config["bowtie2"]["L"],
        D = config["bowtie2"]["D"],
        R = config["bowtie2"]["R"]
    shell: 
        """
        bowtie2 -p {threads} -1 {input.fq1} -2 {input.fq2} -S {output.sam} \
        -x {params.index} -N {params.N} -L {params.L} -D {params.D} -R {params.R} > {log} 2>&1
        """
