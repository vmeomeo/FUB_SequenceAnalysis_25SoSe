# File: rules/map_consensus_patch.smk

rule make_blast_db:
    input:
        fasta = config["candidate_references"]
    output:
        db_files = expand("resources/blastdb/candidates.{ext}", ext=["nhr", "nin", "nsq"])
    conda:
        "../env/blast.yaml"
    shell:
        """
        makeblastdb -in {input.fasta} -dbtype nucl -out resources/blastdb/candidates
        """


rule blast_best_hit:
    input:
        contigs = lambda wc: f"{config['output_dir_path']}/assembly/{wc.sample}/contigs.fasta",
        db = expand("resources/blastdb/candidates.{ext}", ext=["nhr", "nin", "nsq"])
    output:
        best_hit = f"{config['output_dir_path']}/blast/{{sample}}_best_ref.txt"
    conda:
        "../env/blast.yaml"
    shell:
        """
        blastn -query {input.contigs} -db resources/blastdb/candidates -outfmt '6 sseqid bitscore' \
            | sort -k2,2nr | head -n1 | cut -f1 > {output.best_hit}
        """


rule extract_best_reference:
    input:
        best_hit = f"{config['output_dir_path']}/blast/{{sample}}_best_ref.txt",
        references = config["candidate_references"]
    output:
        ref_fasta = f"{config['output_dir_path']}/reference/{{sample}}_best_ref.fasta"
    conda:
        "../env/bioinfo_utils.yaml"
    shell:
        """
        seqtk subseq {input.references} {input.best_hit} > {output.ref_fasta}
        """


rule patch_assembly:
    input:
        contigs = lambda wc: f"{config['output_dir_path']}/polishing/{wc.sample}_polished.fasta" if config.get("enable_polishing", False)
                  else f"{config['output_dir_path']}/assembly/{wc.sample}/contigs.fasta",
        reference = f"{config['output_dir_path']}/reference/{{sample}}_best_ref.fasta"
    output:
        patched = f"{config['output_dir_path']}/patched/{{sample}}_patched.fasta"
    conda:
        "../env/ragtag.yaml"
    shell:
        """
        ragtag.py patch {input.reference} {input.contigs} -o {config['output_dir_path']}/patched/{{sample}}
        cp {config['output_dir_path']}/patched/{{sample}}/ragtag.patch.fasta {output.patched}
        """


rule final_map_with_minimap2:
    input:
        reads = lambda wc: get_preprocessed_reads(wc.sample),
        reference = f"{config['output_dir_path']}/patched/{{sample}}_patched.fasta"
    output:
        sam = f"{config['output_dir_path']}/final_align/sam/{{sample}}.sam"
    threads: workflow.cores
    conda:
        "../env/align_assembly.yaml"
    log:
        f"{config['output_dir_path']}/final_align/logs/align_{{sample}}.log"
    shell:
        """
        minimap2 -t {threads} -ax sr {input.reference} {input.reads[0]} {input.reads[1]} > {output.sam} 2> {log}
        """


rule final_sam_to_sorted_bam:
    input:
        sam = f"{config['output_dir_path']}/final_align/sam/{{sample}}.sam"
    output:
        bam = f"{config['output_dir_path']}/final_align/bam/{{sample}}_sorted.bam"
    threads: workflow.cores
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/final_align/logs/sortbam_{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        """


rule final_index_bam:
    input:
        bam = f"{config['output_dir_path']}/final_align/bam/{{sample}}_sorted.bam"
    output:
        bai = f"{config['output_dir_path']}/final_align/bam/{{sample}}_sorted.bam.bai"
    threads: 1
    conda:
        "../env/samtools.yaml"
    log:
        f"{config['output_dir_path']}/final_align/logs/index_{{sample}}.log"
    shell:
        """
        samtools index {input.bam} > {log} 2>&1
        """


rule final_call_consensus:
    input:
        bam = f"{config['output_dir_path']}/final_align/bam/{{sample}}_sorted.bam",
        bai = f"{config['output_dir_path']}/final_align/bam/{{sample}}_sorted.bam.bai",
        reference = f"{config['output_dir_path']}/patched/{{sample}}_patched.fasta"
    output:
        consensus = f"{config['output_dir_path']}/final_consensus/{{sample}}.fasta"
    conda:
        "../env/samtools.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/final_consensus/logs/consensus_{{sample}}.log"
    shell:
        """
        samtools consensus -aa -f fasta -T {input.reference} -@ {threads} \
            -o {output.consensus} {input.bam} > {log} 2>&1
        """
