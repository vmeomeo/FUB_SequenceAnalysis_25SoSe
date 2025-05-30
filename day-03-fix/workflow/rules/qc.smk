import os

# Test access to a config key
print("Seed mismatches:", config["seed_mismatches"])

# Build a shell command using a config value
command_test = "ls " + config["adapter_file"]
os.system(command_test)

rule fastqc_on_raw:
    input:
        fq = lambda wc: samples.at[wc.sample, f"fq{wc.read_number}"]
    output:
        html = f"{config['output_dir_path']}/qc/fastqc/{{sample}}_{{read_number}}_fastqc.html",
        zip = f"{config['output_dir_path']}/qc/fastqc/{{sample}}_{{read_number}}_fastqc.zip"
    log:
        f"{config['output_dir_path']}/qc/logs/fastqc_on_raw_{{sample}}_{{read_number}}.log"
    threads: 2
    wrapper:
        "v6.2.0/bio/fastqc"



rule trimmomatic:
    input:
        r1 = lambda wc: samples.at[wc.sample, "fq1"],
        r2 = lambda wc: samples.at[wc.sample, "fq2"],
        adapters = config["adapter_file"]
    output:
        r1_paired = f"{config['output_dir_path']}/trimmed/{{sample}}_R1_paired.fastq.gz",
        r1_unpaired = f"{config['output_dir_path']}/trimmed/{{sample}}_R1_unpaired.fastq.gz",
        r2_paired = f"{config['output_dir_path']}/trimmed/{{sample}}_R2_paired.fastq.gz",
        r2_unpaired = f"{config['output_dir_path']}/trimmed/{{sample}}_R2_unpaired.fastq.gz"
    threads: workflow.cores
    conda:
        "../env/trimmomatic.yaml"
    params:
        seed = config["seed_mismatches"],
        pal = config["palindromeClipThreshold"],
        simp = config["simpleClipThreshold"],
        lead = config["leading"],
        trail = config["trailing"],
        sw_lb = config["sliding_window_lb"],
        sw_ub = config["sliding_window_ub"],
        minlen = config["minlen"]
    log:
        f"{config['output_dir_path']}/qc/logs/trimmomatic_{{sample}}.log"
    shell:
        """
        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.r1_paired} {output.r1_unpaired} \
          {output.r2_paired} {output.r2_unpaired} \
          ILLUMINACLIP:{input.adapters}:{params.seed}:{params.pal}:{params.simp} \
          LEADING:{params.lead} TRAILING:{params.trail} \
          SLIDINGWINDOW:{params.sw_lb}:{params.sw_ub} MINLEN:{params.minlen} \
          > {log} 2>&1
        """

rule fastqc_on_processed:
    input:
        "results/trimmed/{sample}_R{read_number}_{paired_or_not}.fastq.gz"
    output:
        html = "results/qc/fastqc/{sample}_R{read_number}_{paired_or_not}_fastqc.html",
        zip = "results/qc/fastqc/{sample}_R{read_number}_{paired_or_not}_fastqc.zip"
    log:
        "results/qc/logs/fastqc_on_processed_{sample}_{read_number}_{paired_or_not}.log"
    threads: 2
    wrapper:
        "v6.2.0/bio/fastqc"

