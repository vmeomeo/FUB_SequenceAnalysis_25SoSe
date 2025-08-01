rule mlst:
    input:
        fasta = "results/assembly/{sample}.fasta"
    output:
        tsv = "results/mlst/{sample}.tsv"
    log:
        "results/mlst/{sample}.log"
    conda:
        "../env/mlst.yaml"
    shell:
        """
        mlst --scheme "Escherichia coli" {input.fasta} > {output.tsv} 2> {log}
        """

rule abricate_card:
    input:
        fasta = "results/assembly/{sample}.fasta"
    output:
        tsv = "results/abricate/card/{sample}.tsv"
    log:
        "results/abricate/card/{sample}.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate --db card {input.fasta} > {output.tsv} 2> {log}
        """

rule abricate_vfdb:
    input:
        fasta = "results/assembly/{sample}.fasta"
    output:
        tsv = "results/abricate/vfdb/{sample}.tsv"
    log:
        "results/abricate/vfdb/{sample}.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate --db vfdb {input.fasta} > {output.tsv} 2> {log}
        """

rule mob_suite:
    input:
        fasta = "results/assembly/{sample}.fasta"
    output:
        directory("results/plasmids/{sample}")
    log:
        "results/plasmids/{sample}/mob_recon.log"
    conda:
        "../env/mob_suite.yaml"
    shell:
        """
        mob_recon --infile {input.fasta} --outdir {output} > {log} 2>&1
        """

rule aggregate_reports:
    input:
        mlst = expand("results/mlst/{sample}.tsv", sample=samples.index),
        card = expand("results/abricate/card/{sample}.tsv", sample=samples.index),
        vfdb = expand("results/abricate/vfdb/{sample}.tsv", sample=samples.index),
        mob  = expand("results/plasmids/{sample}/mobtyper_results.txt", sample=samples.index)
    output:
        summary = "results/summary/combined_report.xlsx"
    conda:
        "../env/summary.yaml"
    log:
        "results/summary/aggregate_reports.log"
    script:
        "workflow/scripts/aggregate_reports.py"

