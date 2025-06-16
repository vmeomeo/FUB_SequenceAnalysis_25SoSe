rule screen_kraken:
    input:
        r1 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_1.fastq",
        r2 = lambda wc: f"{config['output_dir_path']}/cleaned/{wc.sample}_2.fastq"
    output:
        report = f"{config['output_dir_path']}/kraken/{{sample}}.kraken",
        summary = f"{config['output_dir_path']}/kraken/{{sample}}.report"
    conda:
        "../env/kraken_screen.yaml"
    threads: workflow.cores
    params:
        db = config["kraken_db"] 
    log:
        f"{config['output_dir_path']}/kraken/logs/kraken2_{{sample}}.log"
    shell:
        """
        kraken2 --db {params.db} --paired --threads {threads} \
            --report {output.summary} --output {output.report} \
            {input.r1} {input.r2} > {log} 2>&1
        """
