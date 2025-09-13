rule fastp:
    input:
        r1 = lambda wc: samples.at[wc.sample, "fq1"],
        r2 = lambda wc: samples.at[wc.sample, "fq2"]
    output:
        r1_clean = f"{config['output_dir_path']}/cleaned/{{sample}}_r1.fastq.gz",
        r2_clean = f"{config['output_dir_path']}/cleaned/{{sample}}_r2.fastq.gz",
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

rule porechop:
    input:
        lambda wc: samples.at[wc.sample, "fq3"]
    output:
        f"{config['output_dir_path']}/cleaned/{{sample}}_porechop.fastq"
    conda:
        "../env/porechop.yaml"
    log:
        f"{config['output_dir_path']}/qc/logs/porechop_{{sample}}.log"
    shell:
        """
        porechop -i {input} -o {output} > {log} 2>&1
        """

rule nanofilt:
    input:
        f"{config['output_dir_path']}/cleaned/{{sample}}_porechop.fastq"
    output:
        f"{config['output_dir_path']}/filtered/{{sample}}_filtered.fastq"
    conda:
        "../env/nanofilt.yaml"
    log:
        f"{config['output_dir_path']}/qc/logs/nanofilt_{{sample}}.log"
    shell:
        """
        cat {input} | NanoFilt -q 10 -l 1000 > {output} 2> {log}
        """

rule nanoplot:
    input:
        f"{config['output_dir_path']}/filtered/{{sample}}_filtered.fastq"
    output:
        f"{config['output_dir_path']}/qc/nanoplot/{{sample}}/NanoPlot-report.html",
        html_dir = directory(f"{config['output_dir_path']}/qc/nanoplot/{{sample}}")
    conda:
        "../env/nanoplot.yaml"
    log:
        f"{config['output_dir_path']}/qc/logs/nanoplot_{{sample}}.log"
    shell:
        """
        NanoPlot --fastq {input} \
                 --outdir {output.html_dir} \
                 --loglength --N50 > {log} 2>&1
        """

rule multiqc: 
    input:
        fastp_html = expand(f"{config['output_dir_path']}/qc/{{sample}}_fastp.html", sample=samples.index),
        fastp_json = expand(f"{config['output_dir_path']}/qc/{{sample}}_fastp.json", sample=samples.index),
        # nanoplot_summaries = lambda wildcards: (
        #     expand(f"{config['output_dir_path']}/qc/nanoplot/{{sample}}/NanoPlot-report.html", sample=samples.index)
        #     if config.get("long_reads", False) else []
        # )
        nanoplot_summaries = expand(f"{config['output_dir_path']}/qc/nanoplot/{{sample}}/NanoPlot-report.html", sample=samples.index)
    output:
        html = f"{config['output_dir_path']}/qc/multiqc_report.html"
    conda:
        "../env/multiqc.yaml"
    threads: 1
    params:
        outdir = f"{config['output_dir_path']}/qc"
    log:
        f"{config['output_dir_path']}/qc/logs/multiqc.log"
    shell:
        """
        multiqc {input.fastp_json} {input.fastp_html} {input.nanoplot_summaries} -o {params.outdir} > {log} 2>&1
        """


