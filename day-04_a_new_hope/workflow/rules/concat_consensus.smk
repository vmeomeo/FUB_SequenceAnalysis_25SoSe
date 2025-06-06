rule concat_consensus:
    input:
        consensus = expand(
            f"{config['output_dir_path']}/consensus/{{sample}}_{{read_number}}.consensus.fasta",
            sample=samples.index,
            read_number=["1", "2"]
        )
    output:
        combined = f"{config['output_dir_path']}/all_samples/all_samples.fasta"
    conda:
        "../env/concat_consensus.yaml"
    log:
        f"{config['output_dir_path']}/logs/concat_consensus.log"
    shell:
        """
        mkdir -p {config[output_dir_path]}/logs {config[output_dir_path]}/all_samples
        workflow/scripts/select_representative_consensus.py {input.consensus} {output.combined} > {log} 2>&1
        """
