rule bowtie2_align:
    input:
        r1 = lambda wc: get_preprocessed_reads(wc.sample)[0],
        r2 = lambda wc: get_preprocessed_reads(wc.sample)[1],
        contigs = f"{config['output_dir_path']}/assembly/{{sample}}/contigs.fasta"
    output:
        sam = temp(f"{config['output_dir_path']}/polishing/tmp/{{sample}}.sam")
    conda:
        "../env/bowtie2_pre_pol.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/polishing/logs/bowtie2_{{sample}}.log"
    shell:
        """
        bowtie2-build {input.contigs} polishing_contigs_index
        bowtie2 -x polishing_contigs_index \
                -1 {input.r1} -2 {input.r2} -p {threads} \
                -S {output.sam} > {log} 2>&1
        """


rule samtools_process:
    input:
        sam = f"{config['output_dir_path']}/polishing/tmp/{{sample}}.sam"
    output:
        bam = f"{config['output_dir_path']}/polishing/bam/{{sample}}.bam",
        bai = f"{config['output_dir_path']}/polishing/bam/{{sample}}.bam.bai"
    conda:
        "../env/samtools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/polishing/logs/samtools_{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        samtools index {output.bam} > /dev/null 2>> {log}
        """


rule polish_with_pilon:
    input:
        bam = f"{config['output_dir_path']}/polishing/bam/{{sample}}.bam",
        bai = f"{config['output_dir_path']}/polishing/bam/{{sample}}.bam.bai",
        unpolished = f"{config['output_dir_path']}/assembly/{{sample}}/contigs.fasta"
    output:
        polished = f"{config['output_dir_path']}/polishing/{{sample}}_polished.fasta"
    conda:
        "../env/pilon.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/polishing/logs/pilon_{{sample}}.log"
    params:
        outdir = f"{config['output_dir_path']}/polishing"
    shell:
        """
        pilon --genome {input.unpolished} --frags {input.bam} \
              --output {wildcards.sample}_pilon --outdir {params.outdir} \
              --threads {threads} > {log} 2>&1

        mv {params.outdir}/{wildcards.sample}_pilon.fasta {output.polished}
        """
