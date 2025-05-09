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
import os
print(samples.at["sample1", "fq1"])
command_test = "ls " + samples.at["sample1", "fq1"]
os.system(command_test)

rule fastqc_on_raw_1:
    input:
        fq1 = lambda wc: samples.at[wc.sample, "fq1"] 
        # and not "../data/{sample}.fastq.gz" 
    output:
        html = "../results/qc/fastqc/{sample}_1_fastqc.html",
        zip = "../results/qc/fastqc/{sample}_1_fastqc.zip"
    threads: 2
    params:
        outdir=config["output_dir_path"]+"/qc/fastqc"
    wrapper:
        "master/bio/fastqc"

rule fastqc_on_raw_2:
    input:
        fq2 = lambda wc: samples.at[wc.sample, "fq2"] 
        # and not "../data/{sample}.fastq.gz" 
    output:
        html = "../results/qc/fastqc/{sample}_2_fastqc.html",
        zip = "../results/qc/fastqc/{sample}_2_fastqc.zip"
    threads: 2
    params:
        outdir=config["output_dir_path"]+"/qc/fastqc"
    wrapper:
        "master/bio/fastqc"

# rule trimmomatic:
#     input:
#         fq1 = lambda wc: samples.at[wc.sample, "fq1"] 
#         fq2 = lambda wc: samples.at[wc.sample, "fq2"]
#         # and not r1 = "../data/{sample}_1.fastq.gz",
#         #r2 = "../data/{sample}_2.fastq.gz",
#         adapters = config["adapter_file"]
#     output:
#         r1_paired = "../results/trimmed/{sample}_R1_paired.fastq.gz",
#         r1_unpaired = "../results/trimmed/{sample}_R1_unpaired.fastq.gz",
#         r2_paired = "../results/trimmed/{sample}_R2_paired.fastq.gz",
#         r2_unpaired = "../results/trimmed/{sample}_R2_unpaired.fastq.gz"
#     threads: 4
#     conda:
#         "../env/trimmomatic.yaml"
#     shell:
#         """
#         trimmomatic PE -threads config["trimmo_threads"] \
#           {input.r1} {input.r2} \
#           {output.r1_paired} {output.r1_unpaired} \
#           {output.r2_paired} {output.r2_unpaired} \
#           ILLUMINACLIP:{input.adapters}:{config["seed_mismatches"]}:{config["palindromeClipThreshold"]}:{config["simpleClipThreshold"]} \
#           LEADING:{config["leading"]} TRAILING:{config["trailing"]} SLIDINGWINDOW:{config["sliding_window_lb"]}:{config["sliding_window_ub"]} MINLEN:{config["minlen"]}
#         """
# rule fastqc_on_processed:
#     input:
#         "../results/trimmed/{sample}_paired.fastq.gz" 
#     output:
#         html = "qc/fastqc_on_processed/{sample}_fastqc.html",
#         zip = "qc/fastqc_on_processed/{sample}_fastqc.zip"
#     threads: 2
#     wrapper:
#         "0.91.0/bio/fastqc"
