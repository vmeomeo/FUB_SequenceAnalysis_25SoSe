rule featurecounts:
    input:
        annotation = config["annotation"],
        bams = expand(
            f"{config['output_dir_path']}/with_ref_annot/{{sample}}_Aligned.sortedByCoord.out.bam",
            sample=samples.index
        )
    output:
        counts = f"{config['output_dir_path']}/featurecounts/gene_counts.txt",
        summary = f"{config['output_dir_path']}/featurecounts/gene_counts.txt.summary"
    threads: workflow.cores
    conda:
        "../env/subread.yaml"
    log:
        f"{config['output_dir_path']}/featurecounts/logs/featurecounts.log"
    shell:
        """
        featureCounts -T {threads} \
          -a {input.annotation} \
          -o {output.counts} \
          -p -B -C \
          {input.bams} > {log} 2>&1
        """
