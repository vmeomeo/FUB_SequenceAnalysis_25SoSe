# Test
# print(samples.at["ERR024604_tiny", "fq2"]) 

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
        "bowtie2-build -p {threads} -f {input} {config[output_dir_path]}/ref_idx/reference_index "
        "2> {log}"

print(f"Running Bowtie2 with output: {config['output_dir_path']}/sam/{{sample}}.sam")

rule align:
    input: 
        fq1 = "../results/trimmed/{sample}_R1_paired.fastq.gz",
        fq2 = "../results/trimmed/{sample}_R2_paired.fastq.gz",
        index = expand(
            "{output_dir}/ref_idx/reference_index.{ext}",
            output_dir=config["output_dir_path"],
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    output:
        sam=f"{config['output_dir_path']}/sam/{{sample}}.sam"
    threads: workflow.cores
    conda:
        "../env/bowtie2.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_bowtie2_align.log"
    shell: 
        "bowtie2 --no-unal -t --threads {threads} "
        "--np {config[N_penalty]} "
        "--rfg {config[ref_gap_penalty]} " # Reason caused the errors --> solved
        "-x {config[output_dir_path]}/ref_idx/reference_index "
        "-1 {input.fq1} -2 {input.fq2} -S {output.sam}"
        "2> {log}"
