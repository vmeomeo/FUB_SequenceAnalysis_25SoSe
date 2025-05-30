rule qualimap2:
    input:
        bam = f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output:
        report_dir = directory(f"{config['output_dir_path']}/qc/qualimap/{{sample}}_qualimap_report")
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/qc/logs/qualimap_{{sample}}.log"
    conda:
        "../env/qualimap2.yaml"
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {output} -nt {threads} > {log} 2>&1
        """
