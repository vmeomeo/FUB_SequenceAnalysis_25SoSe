rule multiqc:
    input:
        expand("../results/qc/qualimap/{sample}_qualimap_report", sample=samples.index)
    output:
        html = "../results/qc/multiqc_report.html",
        dir = directory("../results/qc/multiqc_data")
    threads: workflow.cores
    conda:
        "../env/multiqc.yaml"
    log:
        "../results/qc/logs/multiqc.log"
    shell:
        """
        multiqc ../results/qc/qualimap -o ../results/qc --threads {threads} > {log} 2>&1
        """

