rule multiqc:
    input:
        qualimap = expand(
            f"{config['output_dir_path']}/qc/qualimap/{{sample}}_qualimap_report",
            sample=samples.index
        ),
        fastqc_raw = expand(
            f"{config['output_dir_path']}/qc/fastqc/{{sample}}_{{read_number}}_fastqc.zip",
            sample=samples.index,
            read_number=["1", "2"]
        ),
        fastqc_trimmed = expand(
            f"{config['output_dir_path']}/qc/fastqc/{{sample}}_R{{read_number}}_{{paired}}_fastqc.zip",
            sample=samples.index,
            read_number=["1", "2"],
            paired=["paired", "unpaired"]
        )
    output:
        html = f"{config['output_dir_path']}/qc/multiqc_report.html",
        dir = directory(f"{config['output_dir_path']}/qc/multiqc_data")
    threads: 1
    conda:
        "../env/multiqc.yaml"
    log:
        f"{config['output_dir_path']}/qc/logs/multiqc.log"
    params:
        outdir = config["output_dir_path"] + "/qc"
    shell:
        """
        multiqc {params.outdir} -o {params.outdir} > {log} 2>&1
        """

