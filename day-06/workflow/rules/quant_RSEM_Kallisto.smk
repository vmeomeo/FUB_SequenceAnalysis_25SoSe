rule kallisto_index_trinity:
    input:
        transcriptome = f"{config['output_dir_path']}/trinity_out/Trinity.fasta"
    output:
        index = f"{config['output_dir_path']}/trinity_out_ka/Trinity.idx"
    log:
        f"{config['output_dir_path']}/trinity_out_ka/kallisto_index_trinity.log"
    conda:
        "../env/kallisto.yaml"
    shell:
        "kallisto index -i {output.index} {input.transcriptome} > {log} 2>&1"

rule kallisto_quant_trinity:
    input:
        index = f"{config['output_dir_path']}/trinity_out_ka/Trinity.idx",
        r1 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz",
        r2 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz"
    output:
        quant = directory(f"{config['output_dir_path']}/quant/trinity_kallisto/{{sample}}")
    log:
        f"{config['output_dir_path']}/quant/trinity_kallisto/{{sample}}/kallisto_quant_trinity.log"
    threads: workflow.cores
    conda:
        "../env/kallisto.yaml"
    shell:
        """
        kallisto quant \
            -i {input.index} \
            -o {output.quant} \
            -b 100 \
            -t {threads} \
            {input.r1} {input.r2} > {log} 2>&1
        """

rule rsem_prepare_trinity_reference:
    input:
        transcriptome = f"{config['output_dir_path']}/trinity_out/Trinity.fasta"
    output:
        f"{config['output_dir_path']}/rsem/Trinity.grp"
    log:
        f"{config['output_dir_path']}/rsem/rsem_prepare_trinity_reference.log"
    conda:
        "../env/rsem.yaml"
    shell:
        """
        rsem-prepare-reference \
            --bowtie2 \
            {input.transcriptome} \
            {config[output_dir_path]}/rsem/Trinity > {log} 2>&1
        """

rule rsem_quant_trinity:
    input:
        r1 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz",
        r2 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz",
        ref = f"{config['output_dir_path']}/rsem/Trinity.grp"
    output:
        genes = f"{config['output_dir_path']}/quant/trinity_rsem/{{sample}}.genes.results",
        isoforms = f"{config['output_dir_path']}/quant/trinity_rsem/{{sample}}.isoforms.results"
    params:
        ref_base = lambda wildcards, input: str(input.ref).removesuffix(".grp"),
        output_prefix = lambda wildcards, output: str(output.genes).removesuffix(".genes.results")
    threads: workflow.cores
    conda:
        "../env/rsem.yaml"
    log:
        f"{config['output_dir_path']}/quant/trinity_rsem/{{sample}}.rsem_quant_trinity.log"
    shell:
        """
        rsem-calculate-expression \
            --paired-end \
            --bowtie2 \
            -p {threads} \
            {input.r1} {input.r2} \
            {params.ref_base} \
            {params.output_prefix} > {log} 2>&1
        """
