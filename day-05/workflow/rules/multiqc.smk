rule multiqc:
    input:
        fastp_reports = expand(f"{config['output_dir_path']}/qc/{{sample}}_{{read_number}}_fastp.json",
                              sample=samples.index, read_number=[1, 2]),
        quast_dirs = expand(f"{config['output_dir_path']}/assembly_qc/{{sample}}_{{read_number}}",
                           sample=samples.index, read_number=[1, 2])
    output:
        report = f"{config['output_dir_path']}/multiqc_report.html"
    conda:
        "../env/multiqc.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/logs/multiqc.log"
    shell:
        """
        multiqc {config['output_dir_path']} -o {config['output_dir_path']} -n multiqc_report.html --threads {threads} > {log} 2>&1
        """
