rule build_contaminant_index:
    input:
        fasta = config["contaminant_fasta"]
    output:
        index_files = expand("resources/contaminant_index.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ])
    conda:
        "../env/bowtie2.yaml"
    threads: 2
    log:
        "resources/logs/build_contaminant_index.log"
    shell:
        """
        bowtie2-build {input.fasta} resources/contaminant_index > {log} 2>&1
        """



rule decontaminate_reads:
    input:
        r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_1.fastq",
        r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_2.fastq",
        index = expand("resources/contaminant_index.{ext}", ext=[
            "1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"
        ])
    output:
        r1_clean = f"{config['output_dir_path']}/decontaminated/{{sample}}_1.fastq",
        r2_clean = f"{config['output_dir_path']}/decontaminated/{{sample}}_2.fastq"
    conda:
        "../env/bowtie2.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/decontaminated/logs/{{sample}}.log"
    shell:
        """
        bowtie2 -x resources/contaminant_index \
            -1 {input.r1} -2 {input.r2} \
            --un-conc results/decontaminated/{wildcards.sample} \
            -p {threads} -S /dev/null 2>> {log}

        mv results/decontaminated/{wildcards.sample}.1 {output.r1_clean}
        mv results/decontaminated/{wildcards.sample}.2 {output.r2_clean}
        """