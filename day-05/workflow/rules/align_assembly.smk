rule align_with_minimap2:
    input:
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}/contigs.fasta",
        reference = "resources/sequence.fasta"
    output:
        f"{config['output_dir_path']}/sam/{{sample}}.sam"
    threads: workflow.cores
    conda:
        "../env/align_assembly.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/align_{{sample}}.log"
    shell:
        """
        minimap2 -t {threads} -ax asm5 {input.reference} {input.assembly} > {output}
        echo "Alignment completed for {wildcards.sample}" > {log}
        """
rule sam_to_sorted_bam:
    input:
        f"{config['output_dir_path']}/sam/{{sample}}.sam"
    output:
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/sam_to_sorted_bam_{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        echo "SAM to sorted BAM done for {wildcards.sample}" > {log}
        """

rule index_bam:
    input: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam.bai"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/alignment/logs/index_bam_{{sample}}.log"
    shell: 
        """
        samtools index {input} {output}
        echo "Indexing done for {wildcards.sample}" > {log}
        """



