# print(lambda wc: samples.at[wc.sample, "fq1"])
print(samples.at["sample1", "fq2"]) 

# rule fastqc_1:
#     input:
#         fq1 = lambda wc: samples.at[wc.sample, "fq1"]
#     output:
#         html1 = config["output_dir_path"] + "/fastqc/{sample}_1_fastqc.html",
#         zip1 = config["output_dir_path"] + "/fastqc/{sample}_1_fastqc.zip"
#     # params:
#     #     outdir = config["output_dir_path"] + "/fastqc"
#     threads: 4
#     conda: "../env/fastqc.yaml"
#     wrapper: "master/bio/fastqc"

# rule fastqc_2:
#     input:
#         fq2 = lambda wc: samples.at[wc.sample, "fq2"]
#     output:
#         html2 = config["output_dir_path"] + "/fastqc/{sample}_2_fastqc.html",
#         zip2 = config["output_dir_path"] + "/fastqc/{sample}_2_fastqc.zip"
#     # params:
#     #     outdir = config["output_dir_path"] + "/fastqc"
#     threads: 4
#     conda: "../env/fastqc.yaml"
#     wrapper: "master/bio/fastqc"

rule fastqc:
    input:
        fq1 = lambda wc: samples.at[wc.sample, "fq1"],
        fq2 = lambda wc: samples.at[wc.sample, "fq2"]
    output:
        html1 = config["output_dir_path"] + "/fastqc/{sample}_1_fastqc.html",
        html2 = config["output_dir_path"] + "/fastqc/{sample}_2_fastqc.html",
        zip1 = config["output_dir_path"] + "/fastqc/{sample}_1_fastqc.zip",
        zip2 = config["output_dir_path"] + "/fastqc/{sample}_2_fastqc.zip"
    params:
        outdir = config["output_dir_path"] + "/fastqc"
    threads: 4
    conda: "../env/fastqc.yaml"
    shell: 
        "fastqc -o {params.outdir} {input.fq1} {input.fq2}"