rule align_with_minimap2:
    input:
        assembly = lambda wc: f"{config['output_dir_path']}/polishing/{wc.sample}_polished.fasta"
        if config.get("enable_polishing", False)
        else f"{config['output_dir_path']}/assembly/{wc.sample}/contigs.fasta",
        reference = config["reference_genome"]
    output:
        sam = f"{config['output_dir_path']}/alignment/sam/{{sample}}.sam"
    threads: workflow.cores
    conda:
        "../env/align_assembly.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/align_{{sample}}.log"
    shell:
        """
        minimap2 -t {threads} -ax asm5 {input.reference} {input.assembly} > {output.sam} 2> {log}
        """


rule sam_to_sorted_bam:
    input:
        sam = f"{config['output_dir_path']}/alignment/sam/{{sample}}.sam"
    output:
        bam = f"{config['output_dir_path']}/alignment/bam/{{sample}}_sorted.bam"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/sortbam_{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        """


rule index_bam:
    input:
        bam = f"{config['output_dir_path']}/alignment/bam/{{sample}}_sorted.bam"
    output:
        bai = f"{config['output_dir_path']}/alignment/bam/{{sample}}_sorted.bam.bai"
    threads: 1
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/index_{{sample}}.log"
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """


rule call_consensus:
    input:
        bam = f"{config['output_dir_path']}/alignment/bam/{{sample}}_sorted.bam",
        bai = f"{config['output_dir_path']}/alignment/bam/{{sample}}_sorted.bam.bai",
        reference = config["reference_genome"]
    output:
        consensus = f"{config['output_dir_path']}/consensus/{{sample}}.consensus.fasta"
    conda:
        "../env/samtools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/consensus/logs/consensus_{{sample}}.log"
    shell:
        """
        samtools consensus -aa -f fasta -T {input.reference} -@ {threads} \
            -o {output.consensus} {input.bam} > {log} 2>&1
        """


rule concat_consensus:
    input:
        expand(
            f"{config['output_dir_path']}/consensus/{{sample}}.consensus.fasta",
            sample=samples.index
        )
    output:
        combined = f"{config['output_dir_path']}/all_samples/all_samples.fasta"
    conda:
        "../env/concat_consensus.yaml"
    log:
        f"{config['output_dir_path']}/logs/concat_consensus.log"
    shell:
        """
        workflow/scripts/select_representative_consensus.py {input} {output.combined} > {log} 2>&1
        """