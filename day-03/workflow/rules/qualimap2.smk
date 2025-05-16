rule qualimap2:
    input:
        bam=f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam"
    output:
        report_dir=directory("../results/qc/qualimap/{sample}_qualimap_report")
    threads: workflow.cores
    log:
        "../results/qc/logs/qualimap_{sample}.log"
    conda:
        "../env/qualimap2.yaml"
    params:
        nt=workflow.cores
    wrapper:
        "master/bio/qualimap/bamqc"
    # shell:
    #     """
    #     qualimap bamqc -bam {input.bam} -outdir {output.report_dir} --num-threads {threads} > {log} 2>&1
    #     """
