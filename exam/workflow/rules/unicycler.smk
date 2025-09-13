def get_short_read1(wildcards):
    if config.get("contaminant_fasta", ""):
        return f"{config['output_dir_path']}/decontaminated/{wildcards.sample}_r1.fastq"
    return f"{config['output_dir_path']}/cleaned/{wildcards.sample}_r1.fastq.gz"

def get_short_read2(wildcards):
    if config.get("contaminant_fasta", ""):
        return f"{config['output_dir_path']}/decontaminated/{wildcards.sample}_r2.fastq"
    return f"{config['output_dir_path']}/cleaned/{wildcards.sample}_r2.fastq.gz"

def get_long_read(wildcards):
    if config.get("long_reads", False):
        return f"{config['output_dir_path']}/filtered/{wildcards.sample}_filtered.fastq"
    return None  # Explicitly return None if no long reads
rule unicycler:
    input:
        r1 = get_short_read1,
        r2 = get_short_read2,
        l = get_long_read,
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
            -l {input.l} \
            -o {params.outdir}/assembly/{wildcards.sample}/temp -t {threads} \
            &> {log}
        mv {params.outdir}/assembly/{wildcards.sample}/temp/assembly.fasta {output}
        """