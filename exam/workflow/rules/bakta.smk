rule download_bakta_db:
    output:
        directory(config["bakta_db_path"])
    log:
        f"{config['output_dir_path']}/logs/bakta/download_bakta_db.log"
    conda:
        "../env/bakta.yaml"
    shell:
        "bakta_db download --output {output} --type light &> {log}"

rule bakta:
    input:
        f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta"
    output:
        gff = f"{config['output_dir_path']}/annotation/{{sample}}/{{sample}}.gff3",
        faa = f"{config['output_dir_path']}/annotation/{{sample}}/{{sample}}.faa",
        ffn = f"{config['output_dir_path']}/annotation/{{sample}}/{{sample}}.ffn"
    log:
        f"{config['output_dir_path']}/logs/bakta/{{sample}}.log"
    threads: workflow.cores
    conda:
        "../env/bakta.yaml"
    params:
        db = config['bakta_db_path'],
        outdir = config['output_dir_path']
    shell:
        "bakta --db {params.db} {input} "
        "--output {params.outdir}/annotation/{wildcards.sample} "
        "--threads {threads} --force &> \"{log}\""


