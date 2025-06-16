rule screen:
    input:
        expand(f"{config['output_dir_path']}/kraken/{{sample}}.report", sample=samples.index),
        f"{config['output_dir_path']}/kraken/multiqc_report.html"
