rule call_consensus:
    input:
        bam = f"{config['output_dir_path']}/bam_sorted/{{sample}}_{{read_number}}_sorted.bam",
        bai = f"{config['output_dir_path']}/bam_sorted/{{sample}}_{{read_number}}_sorted.bam.bai",
        reference = "resources/sequence.fasta"
    output:
        vcf = f"{config['output_dir_path']}/consensus/{{sample}}_{{read_number}}.vcf.gz",
        vcf_index = f"{config['output_dir_path']}/consensus/{{sample}}_{{read_number}}.vcf.gz.tbi",
        consensus = f"{config['output_dir_path']}/consensus/{{sample}}_{{read_number}}.consensus.fasta"
    conda:
        "../env/bcftools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/consensus/logs/consensus_{{sample}}_{{read_number}}.log"
    shell:
        """
        mkdir -p {config[output_dir_path]}/consensus/logs

        bcftools mpileup -f {input.reference} {input.bam} |
        bcftools call -c -Oz -o {output.vcf}

        tabix -p vcf {output.vcf}

        bcftools consensus -f {input.reference} {output.vcf} > {output.consensus}

        echo "Consensus sequence generated for {wildcards.sample}_{wildcards.read_number}" > {log}
        """
