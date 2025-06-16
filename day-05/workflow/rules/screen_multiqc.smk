rule screen_multiqc:
    input:
        reports = expand(f"{config['output_dir_path']}/kraken/{{sample}}.report", sample=samples.index)
    output:
        html = f"{config['output_dir_path']}/kraken/multiqc_report.html",
        dir = directory(f"{config['output_dir_path']}/kraken/multiqc_data")
    conda:
        "../env/multiqc_screen.yaml"
    threads: 1
    log:
        f"{config['output_dir_path']}/kraken/logs/multiqc_kraken.log"
    params:
        outdir = f"{config['output_dir_path']}/kraken"
    shell:
        """
        multiqc {input.reports} -o {params.outdir} > {log} 2>&1
        """