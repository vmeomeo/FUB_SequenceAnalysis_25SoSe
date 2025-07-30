rule merge_all_reads:
    input:
        r1 = samples.index.map(lambda s: f"{config['output_dir_path']}/cleaned/{s}_r1.fastq.gz").tolist(),
        r2 = samples.index.map(lambda s: f"{config['output_dir_path']}/cleaned/{s}_r2.fastq.gz").tolist()
    output:
        r1 = f"{config['output_dir_path']}/merged/all_r1.fastq.gz",
        r2 = f"{config['output_dir_path']}/merged/all_r2.fastq.gz"
    log:
        log_r1 = f"{config['output_dir_path']}/logs/merge_all_r1.log",
        log_r2 = f"{config['output_dir_path']}/logs/merge_all_r2.log"
    threads: workflow.cores
    shell:
        """        
        echo "Merging R1 files: {input.r1}" > {log.log_r1}
        zcat {input.r1} | gzip > {output.r1} 2>> {log.log_r1} || exit 1

        echo "Merging R2 files: {input.r2}" > {log.log_r2}
        zcat {input.r2} | gzip > {output.r2} 2>> {log.log_r2} || exit 1
        """

rule trinity_merged:
    input:
        left = f"{config['output_dir_path']}/merged/all_r1.fastq.gz",
        right = f"{config['output_dir_path']}/merged/all_r2.fastq.gz"
    output:
        outdir = directory(f"{config['output_dir_path']}/trinity_out"),
        fasta = f"{config['output_dir_path']}/trinity_out/Trinity.fasta"
    threads: workflow.cores
    conda:
        "../env/trinity.yaml"
    log:
        f"{config['output_dir_path']}/trinity_out/trinity.log"
    shell:
        """
        Trinity \
            --seqType fq \
            --max_memory 50G \
            --left {input.left} \
            --right {input.right} \
            --CPU {threads} \
            --output {output.outdir} \
            --normalize_reads \
            > {log} 2>&1
        """



