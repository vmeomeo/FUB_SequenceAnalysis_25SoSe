rule sam_to_bam:
    input: 
        f"{config['output_dir_path']}/sam/{{sample}}.sam"
    output: 
        f"{config['output_dir_path']}/bam/{{sample}}.bam"
    threads: config["threads"]
    conda:
        "../env/samtools.yaml"
    shell: "samtools view -bS {input} > {output}"

rule sort_bam:
    input: 
        f"{config['output_dir_path']}/bam/{{sample}}.bam"
    output: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    threads: config["threads"]
    conda:
        "../env/samtools.yaml"
    shell: "samtools sort {input} > {output}"

rule index_bam:
    input: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam.bai"
    threads: config["threads"]
    conda:
        "../env/samtools.yaml"
    shell: "samtools index {input} {output}"

rule stats_bam:
    input: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output: 
        f"{config['output_dir_path']}/stats/{{sample}}.stats"
    threads: config["threads"]
    conda:
        "../env/samtools.yaml"
    shell: "samtools idxstats {input} > {output}"