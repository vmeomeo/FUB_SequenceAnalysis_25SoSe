# rule ragtag:
#     input:
#         assembly = f"{config['output_dir_path']}/assembly/{{sample}}_{{read_number}}/contigs.fasta",
#         reference = config["reference_genome"]
#     output:
#         scaffolded = f"{config['output_dir_path']}/scaffolded/{{sample}}_{{read_number}}.scaffolded.fasta"
#     log:
#         f"{config['output_dir_path']}/logs/ragtag_scaffold_{{sample}}_{{read_number}}.log"
#     conda:
#         "../env/ragtag.yaml"
#     shell:
#         """
#         mkdir -p $(dirname {output.scaffolded})
#         ragtag.py scaffold -u {input.reference} {input.assembly} -o $(dirname {output.scaffolded}) > {log} 2>&1
#         cp $(dirname {output.scaffolded})/ragtag.scaffold.fasta {output.scaffolded}
#         """
rule ragtag:
    input:
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}/final.contigs.fa",
        reference = config["reference_genome"]
    output:
        scaffolded = f"{config['output_dir_path']}/scaffolded/{{sample}}.scaffolded.fasta"
    log:
        f"{config['output_dir_path']}/scaffolded/logs/ragtag_scaffold_{{sample}}.log"
    conda:
        "../env/ragtag.yaml"
    shell:
        """
        mkdir -p $(dirname {output.scaffolded})
        ragtag.py scaffold -u {input.reference} {input.assembly} -o $(dirname {output.scaffolded}) > {log} 2>&1
        cp $(dirname {output.scaffolded})/ragtag.scaffold.fasta {output.scaffolded}
        """
