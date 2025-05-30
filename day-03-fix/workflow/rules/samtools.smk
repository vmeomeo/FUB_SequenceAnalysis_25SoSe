rule sam_to_bam:
    input: 
        f"{config['output_dir_path']}/sam/{{sample}}.sam"
    output: 
        f"{config['output_dir_path']}/bam/{{sample}}.bam"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_sam_to_bam.log"
    shell: 
        "samtools view -@ {threads} -bS {input} > {output} 2> {log}"


rule sort_bam:
    input: 
        f"{config['output_dir_path']}/bam/{{sample}}.bam"
    output: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_sort_bam.log"
    shell: 
        "samtools sort -@ {threads} {input} -o {output} 2> {log}"


rule index_bam:
    input: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output: 
        f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam.bai"
    threads: 1
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_index_bam.log"
    shell: 
        "samtools index {input} {output} 2> {log}"
#this above rule index_bam is only able to use 1 thread

rule stats_bam:
    input: 
        bam=f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam",
        index=f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam.bai"
    output: 
        f"{config['output_dir_path']}/stats/{{sample}}.stats"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/logs/{{sample}}_stats_bam.log"
    shell: 
        "samtools idxstats {input.bam} > {output} 2> {log}"
