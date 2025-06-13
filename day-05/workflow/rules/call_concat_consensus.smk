rule call_consensus:
    input:
        bam = f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam",
        bai = f"{config['output_dir_path']}/bam_sorted/{{sample}}_sorted.bam.bai",
        reference = config["reference_genome"]
    output:
        consensus = f"{config['output_dir_path']}/consensus/{{sample}}.consensus.fasta"
    conda:
        "../env/samtools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/consensus/logs/consensus_{{sample}}.log"
    shell:
        """
        samtools consensus -aa -f fasta -T {input.reference} -@ {threads} -o {output.consensus} {input.bam}
        echo "Consensus sequence generated for {wildcards.sample}" > {log}
        """
        
rule concat_consensus:
    input:
        consensus = expand(
            f"{config['output_dir_path']}/consensus/{{sample}}.consensus.fasta",
            sample=samples.index
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
