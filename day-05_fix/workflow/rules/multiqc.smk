rule multiqc:
    input:
        fastp_json = expand(
            f"{config['output_dir_path']}/qc/{{sample}}_fastp.json",
            sample=samples.index
        ),
        fastp_html = expand(
            f"{config['output_dir_path']}/qc/{{sample}}_fastp.html",
            sample=samples.index
        ),
        quast_dirs = expand(
            f"{config['output_dir_path']}/assembly_qc/{{sample}}",
            sample=samples.index
        )
    output:
        html = f"{config['output_dir_path']}/qc/multiqc_report.html"
    conda:
        "../env/multiqc.yaml"
    threads: 1
    params:
        outdir = f"{config['output_dir_path']}/qc"
    log:
        f"{config['output_dir_path']}/logs/multiqc.log"
    shell:
        """
        multiqc {input.fastp_json} {input.fastp_html} {input.quast_dirs} -o {params.outdir} > {log} 2>&1
        """