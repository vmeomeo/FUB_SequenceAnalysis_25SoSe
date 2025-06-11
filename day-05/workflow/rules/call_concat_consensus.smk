rule call_consensus:
    input:
        bam = f"{config['output_dir_path']}/bam_sorted/{{sample}}_{{read_number}}_sorted.bam",
        bai = f"{config['output_dir_path']}/bam_sorted/{{sample}}_{{read_number}}_sorted.bam.bai",
        reference = config["reference_genome"]
    output:
        consensus = f"{config['output_dir_path']}/consensus/{{sample}}_{{read_number}}.consensus.fasta"
    conda:
        "../env/samtools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/consensus/logs/consensus_{{sample}}_{{read_number}}.log"
    shell:
        """
        samtools mpileup {config[samtools_mpileup_flags]} -@ {threads} -f {input.reference} {input.bam} | \
        samtools consensus {config[samtools_consensus_flags]} -@ {threads} -o {output.consensus} -

        echo "Consensus sequence generated for {wildcards.sample}_{wildcards.read_number}" > {log}
        """
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
        workflow/scripts/select_representative_consensus.py {input.consensus} {output.combined} > {log} 2>&1
        """
