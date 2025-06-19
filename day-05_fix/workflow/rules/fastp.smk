print(expand(samples.at["sample_22", "fq1"]))

rule fastp:
    input:
        r1 = lambda wc: samples.at[wc.sample, "fq1"],
        r2 = lambda wc: samples.at[wc.sample, "fq2"]
    output:
        r1_clean = f"{config['output_dir_path']}/cleaned/{{sample}}_1.fastq",
        r2_clean = f"{config['output_dir_path']}/cleaned/{{sample}}_2.fastq",
        html = f"{config['output_dir_path']}/qc/{{sample}}_fastp.html",
        json = f"{config['output_dir_path']}/qc/{{sample}}_fastp.json"
    params:
        flags = config['fastp_flags']
    conda:
        "../env/fastp.yaml"
    threads: 4
    log:
        f"{config['output_dir_path']}/qc/logs/fastp_{{sample}}.log"
    shell:
        """
        fastp {params.flags} -i {input.r1} -I {input.r2} \
              -o {output.r1_clean} -O {output.r2_clean} \
              -j {output.json} -h {output.html} > {log} 2>&1
        """
