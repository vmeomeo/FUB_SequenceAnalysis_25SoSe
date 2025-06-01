rule concat:
    input:
        scaffolds = expand(
            f"{config['output_dir_path']}/scaffolded/{{sample}}_{{read_number}}.scaffolded.fasta",
            sample=samples.index,
            read_number=["1", "2"]
        )
    output:
        combined = f"{config['output_dir_path']}/all_samples/all_samples.fasta"
    conda:
        "../env/concat.yaml"
    log:
        f"{config['output_dir_path']}/logs/collect_representative_sequences.log"
    shell:
        """
        workflow/scripts/select_representative.py {input.scaffolds} {output.combined}
        """