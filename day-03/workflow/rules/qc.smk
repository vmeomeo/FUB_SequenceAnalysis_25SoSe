# rule fastqc:
#     input:
#         fq1 = lambda wc: samples.at[wc.sample, "fq1"],
#         fq2 = lambda wc: samples.at[wc.sample, "fq2"]
#     output: 
#     params:
#     threads: 8
#     conda: "../env/fastqc.yaml"
#     wrapper: 
#         "fastqc {input.fq1} {input.fq2}"

#         fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
#            [-c contaminant file] seqfile1 .. seqfileN
# import os
# print(samples.at["sample1", "fq1"])
# command_test = "ls " + samples.at["sample1", "fq1"]
# os.system(command_test)
import os

# Test access to a config key
print("Seed mismatches:", config["seed_mismatches"])

# Build a shell command using a config value
command_test = "ls " + config["adapter_file"]
os.system(command_test)

# rule fastqc_on_raw_1:
#     input:
#         fq1 = lambda wc: samples.at[wc.sample, "fq1"] 
#         # and not "../data/{sample}.fastq.gz" 
#     output:
#         html = "../results/qc/fastqc/{sample}_1_fastqc.html",
#         zip = "../results/qc/fastqc/{sample}_1_fastqc.zip"
#     threads: 2
#     params:
#         outdir=config["output_dir_path"]+"/qc/fastqc"
#     wrapper:
#         "master/bio/fastqc"

rule fastqc_on_raw:
    input:
        fq = lambda wc: samples.at[wc.sample, f"fq{wc.read_number}"]
    output:
        html = "../results/qc/fastqc/{sample}_{read_number}_fastqc.html",
        zip = "../results/qc/fastqc/{sample}_{read_number}_fastqc.zip"
    threads: workflow.cores
    params:
        outdir = config["output_dir_path"] + "/qc/fastqc"
    log:
        "../results/qc/logs/fastqc_on_raw_{sample}_{read_number}.log"
    conda:
        "../env/fastqc.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.fq} > {log} 2>&1
        """


rule trimmomatic:
    input:
        r1 = lambda wc: samples.at[wc.sample, "fq1"],
        r2 = lambda wc: samples.at[wc.sample, "fq2"],
        adapters = config["adapter_file"]
    output:
        r1_paired = "../results/trimmed/{sample}_R1_paired.fastq.gz",
        r1_unpaired = "../results/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_paired = "../results/trimmed/{sample}_R2_paired.fastq.gz",
        r2_unpaired = "../results/trimmed/{sample}_R2_unpaired.fastq.gz"
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
        "../results/qc/logs/trimmomatic_{sample}.log"
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
        i1= "../results/trimmed/{sample}_R{read_number}_{paired_or_not}.fastq.gz"
    output:
        html = "../results/qc/fastqc/{sample}_{read_number}_fastqc_on_processed_{paired_or_not}.html",
        zip = "../results/qc/fastqc/{sample}_{read_number}_fastqc_on_processed_{paired_or_not}.zip"
    threads: workflow.cores
    params:
        outdir = config["output_dir_path"] + "/qc/fastqc"
    log:
        "../results/qc/logs/fastqc_on_processed_{sample}_{read_number}_{paired_or_not}.log"
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.i1} > {log} 2>&1
        """


