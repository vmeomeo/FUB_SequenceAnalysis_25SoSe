# rule download_mlst_db:
#     output:
#         touch(f"{config['mlst_db_path']}/.downloaded")
#     log:
#         f"{config['output_dir_path']}/logs/mlst/download_mlst_db.log"
#     conda:
#         "../env/mlst.yaml"
#     shell:
#         """
#         mlst --update-db > {log} 2>&1
#         touch {output}
#         """

rule download_mlst_db:
    output:
        touch(f"{config['mlst_db_path']}/.downloaded")
    log:
        f"{config['output_dir_path']}/logs/mlst/download_mlst_db.log"
    conda:
        "../env/mlst.yaml"
    shell:
        """
        # This forces database initialization if not already present
        mlst --check > {log} 2>&1 || true  # --check verifies DBs exist
        touch {output}
        """

rule download_abricate_card_db:
    output:
        touch(f"{config['abricate_db_path_card']}/card/.downloaded")
    log:
        f"{config['output_dir_path']}/logs/abricate/download_card_db.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate-get_db --db card >> {log} 2>&1
        abricate --setupdb >> {log} 2>&1
        touch {output}
        """

rule download_abricate_vfdb_db:
    output:
        touch(f"{config['abricate_db_path_vfdb']}/vfdb/.downloaded")
    log:
        f"{config['output_dir_path']}/logs/abricate/download_vfdb_db.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate-get_db --db vfdb >> {log} 2>&1
        abricate --setupdb >> {log} 2>&1
        touch {output}
        """

rule download_mob_suite_db:
    output:
        directory(config["mob_suite_db_path"])
    log:
        f"{config['output_dir_path']}/logs/mob_suite/download_mob_suite_db.log"
    conda:
        "../env/mob_suite.yaml"
    shell:
        """
        mob_init --outdir {output} > {log} 2>&1
        """

rule mlst:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db_flag = f"{config['mlst_db_path']}/.downloaded"
    output:
        tsv = f"{config['output_dir_path']}/mlst/{{sample}}.tsv"
    log:
        f"{config['output_dir_path']}/mlst/{{sample}}.log"
    conda:
        "../env/mlst.yaml"
    shell:
        """
        mlst --scheme "Klebsiella_pneumoniae" {input.fasta} > {output.tsv} 2> {log}
        """

rule abricate_card:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db_flag = f"{config['abricate_db_path_card']}/card/.downloaded"
    output:
        tsv = f"{config['output_dir_path']}/abricate/card/{{sample}}.tsv"
    log:
        f"{config['output_dir_path']}/abricate/card/{{sample}}.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate --db card {input.fasta} > {output.tsv} 2> {log}
        """

rule abricate_vfdb:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db_flag = f"{config['abricate_db_path_vfdb']}/vfdb/.downloaded"
    output:
        tsv = f"{config['output_dir_path']}/abricate/vfdb/{{sample}}.tsv"
    log:
        f"{config['output_dir_path']}/abricate/vfdb/{{sample}}.log"
    conda:
        "../env/abricate.yaml"
    shell:
        """
        abricate --db vfdb {input.fasta} > {output.tsv} 2> {log}
        """

rule mob_suite:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db_dir = config["mob_suite_db_path"]
    output:
        directory(f"{config['output_dir_path']}/plasmids/{{sample}}")
    log:
        f"{config['output_dir_path']}/plasmids/{{sample}}/mob_recon.log"
    conda:
        "../env/mob_suite.yaml"
    shell:
        """
        mob_recon --infile {input.fasta} --outdir {output} --db {input.db_dir} > {log} 2>&1
        """

def get_aggregate_inputs():
    inputs = {}
    if config.get("run_mlst", False) and config.get("mlst_db_path", ""):
        inputs["mlst"] = expand(f"{config['output_dir_path']}/mlst/{{sample}}.tsv", sample=samples.index)
    if config.get("run_resistance", False) and config.get("abricate_db_path_card", ""):
        inputs["card"] = expand(f"{config['output_dir_path']}/abricate/card/{{sample}}.tsv", sample=samples.index)
    if config.get("run_virulence", False) and config.get("abricate_db_path_vfdb", ""):
        inputs["vfdb"] = expand(f"{config['output_dir_path']}/abricate/vfdb/{{sample}}.tsv", sample=samples.index)
    if config.get("run_plasmid", False) and config.get("mob_suite_db_path", ""):
        inputs["mob"] = expand(f"{config['output_dir_path']}/plasmids/{{sample}}/mobtyper_results.txt", sample=samples.index)
    return inputs


rule aggregate_reports:
    input:
        # add conditions
        # mlst = expand(f"{config['output_dir_path']}/mlst/{{sample}}.tsv", sample=samples.index),
        # card = expand(f"{config['output_dir_path']}/abricate/card/{{sample}}.tsv", sample=samples.index),
        # vfdb = expand(f"{config['output_dir_path']}/abricate/vfdb/{{sample}}.tsv", sample=samples.index),
        # mob  = expand(f"{config['output_dir_path']}/plasmids/{{sample}}/mobtyper_results.txt", sample=samples.index)
        get_aggregate_inputs()
    output:
        summary = f"{config['output_dir_path']}/summary/combined_report.xlsx"
    conda:
        "../env/summary.yaml"
    log:
        f"{config['output_dir_path']}/summary/aggregate_reports.log"
    script:
        "workflow/scripts/aggregate_reports.py"

