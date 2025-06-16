rule build_contaminant_index:
    input:
        fasta = config.get("contaminant_fasta", "")
    output:
        touch("resources/.contaminant_index_placeholder"),  # Creates empty file
        # Original index files as optional_outputs
        *expand("resources/contaminant_index.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    conda:
        "../env/bowtie2.yaml"
    threads: 2
    log:
        "resources/logs/build_contaminant_index.log"
    shell:
        """
        if [ -z "{input.fasta}" ]; then
            echo "No contaminant_fasta provided - skipping index creation" > {log}
            exit 0
        fi
        bowtie2-build {input.fasta} resources/contaminant_index > {log} 2>&1
        """

rule decontaminate_reads:
    input:
        r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_1.fastq",
        r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_2.fastq",
        index_placeholder = "resources/.contaminant_index_placeholder"
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
        if [ ! -f "resources/contaminant_index.1.bt2" ]; then
            echo "No contaminant index found - copying cleaned reads as decontaminated" > {log}
            cp {input.r1} {output.r1_clean}
            cp {input.r2} {output.r2_clean}
        else
            bowtie2 -x resources/contaminant_index -1 {input.r1} -2 {input.r2} \
                --un-conc {config[output_dir_path]}/decontaminated/{{sample}}.fastq \
                -p {threads} -S /dev/null > {log} 2>&1
            mv {config[output_dir_path]}/decontaminated/{{sample}}.1.fastq {output.r1_clean}
            mv {config[output_dir_path]}/decontaminated/{{sample}}.2.fastq {output.r2_clean}
        fi
        """