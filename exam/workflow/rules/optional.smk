with open("config/mlst_default_db.txt","r") as f:
    string = f.read()
mlst_db = string.split(" ")
if config['mlst_scheme'] not in mlst_db:
    raise ValueError(f"MLST scheme '{config['mlst_scheme']}' not found in the default database. Please provide a valid scheme or find the data for your organism on https://pubmlst.org/. You can use the 'mlst_db_path' in the config to point to a custom database. If not provided, the default database will be used. If no data is available for your organism, you can set 'run_mlst' to false in the config file to skip this step.")


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

# rule download_mob_suite_db:
#     output:
#         directory(config["mob_suite_db_path"])
#     log:
#         f"{config['output_dir_path']}/logs/mob_suite/download_mob_suite_db.log"
#     conda:
#         "../env/mob_suite.yaml"
#     shell:
#         """
#         mob_init --outdir {output} > {log} 2>&1
#         """

rule download_plasmidfinder_db:
    output:
        directory(config["plasmidfinder_db_path"])
    log:
        f"{config['output_dir_path']}/logs/plasmidfinder/download_db.log"
    conda:
        "../env/plasmidfinder.yaml"
    shell:
        """
        plasmidfinder.py --download -o {output} > {log} 2>&1
        """



# rule mlst:
#     input:
#         fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
#         db_flag = f"{config['mlst_db_path']}/.downloaded"
#     output:
#         tsv = f"{config['output_dir_path']}/mlst/{{sample}}.tsv"
#     log:
#         f"{config['output_dir_path']}/mlst/{{sample}}.log"
#     conda:
#         "../env/mlst.yaml"
#     params:
#         scheme = config['mlst_scheme']
#     shell:
#         """
#         mlst --scheme "{params.scheme}" {input.fasta} > {output.tsv} 2> {log}
#         """

rule mlst:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta"
    output:
        tsv = f"{config['output_dir_path']}/mlst/{{sample}}.tsv"
    log:
        f"{config['output_dir_path']}/mlst/{{sample}}.log"
    conda:
        "../env/mlst.yaml"
    params:
        scheme = config['mlst_scheme'],
        db_dir = config.get('mlst_db_path', '') 
    shell:
        """
        mlst --scheme "{params.scheme}" \
             {input.fasta} > {output.tsv} 2> {log}
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

# rule mob_suite:
#     input:
#         fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
#         db_dir = config["mob_suite_db_path"]
#     output:
#         directory(f"{config['output_dir_path']}/plasmids/{{sample}}")
#     log:
#         f"{config['output_dir_path']}/plasmids/{{sample}}/mob_recon.log"
#     conda:
#         "../env/mob_suite.yaml"
#     shell:
#         """
#         mob_recon -u --infile {input.fasta} --outdir {output} --db {input.db_dir} > {log} 2>&1
#         """

rule plasmidfinder:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db = config["plasmidfinder_db_path"]
    output:
        result_dir = directory(f"{config['output_dir_path']}/plasmidfinder/{{sample}}")
    log:
        f"{config['output_dir_path']}/plasmidfinder/{{sample}}/plasmidfinder.log"
    conda:
        "../env/plasmidfinder.yaml"
    shell:
        """
        plasmidfinder.py -i {input.fasta} -o {output.result_dir} -p {input.db} > {log} 2>&1
        """



def get_aggregate_inputs():
    inputs = []
    
    # MLST results
    if config.get("run_mlst", False) and config.get("mlst_db_path", ""):
        inputs.extend(expand(f"{config['output_dir_path']}/mlst/{{sample}}.tsv", sample=get_included_samples()))
    
    # CARD resistance genes
    if config.get("run_resistance", False) and config.get("abricate_db_path_card", ""):
        inputs.extend(expand(f"{config['output_dir_path']}/abricate/card/{{sample}}.tsv", sample=get_included_samples()))
    
    # VFDB virulence factors
    if config.get("run_virulence", False) and config.get("abricate_db_path_vfdb", ""):
        inputs.extend(expand(f"{config['output_dir_path']}/abricate/vfdb/{{sample}}.tsv", sample=get_included_samples()))
    
    # PlasmidFinder results
    if config.get("run_plasmidfinder", False) and config.get("plasmidfinder_db_path", ""):
        inputs.extend(directory(expand(f"{config['output_dir_path']}/plasmidfinder/{{sample}}", sample=get_included_samples())))
    
    return inputs

rule aggregate_reports:
    input:
        get_aggregate_inputs()
    output:
        summary = f"{config['output_dir_path']}/summary/combined_report.xlsx"
    conda:
        "../env/summary.yaml"
    log:
        f"{config['output_dir_path']}/summary/aggregate_reports.log"
    script:
        "workflow/scripts/aggregate_reports.py"