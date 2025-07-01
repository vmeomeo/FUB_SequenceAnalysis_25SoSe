rule decompress_fasta:
    input:
        config['ref_genome']
    output:
        config['ref_genome'].replace('.gz', '')
    shell:
        "gunzip -c {input} > {output}"

rule decompress_gtf:
    input:
        config['annotation']
    output:
        config['annotation'].replace('.gz', '')
    shell:
        "gunzip -c {input} > {output}"

rule genome_indexing:
    input:
        fa = config['ref_genome'].replace('.gz', ''),
        gtf = config['annotation'].replace('.gz', '')
    output:
        directory("resources/index")
    threads: 4
    conda:
        "../env/star.yaml"
    log:
        f"{config['output_dir_path']}/indexing/logs/STAR_genome_indexing.log"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99 \
        > {log} 2>&1
        """
rule genome_mapping:
    input:
        idx = "resources/index",
        r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_r1.fastq.gz",
        r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_r2.fastq.gz"
    output:
        bam = f"{config['output_dir_path']}/with_ref_annot/{{sample}}_Aligned.sortedByCoord.out.bam",
        counts = f"{config['output_dir_path']}/with_ref_annot/{{sample}}_ReadsPerGene.out.tab"
    threads: 4
    conda:
        "../env/star.yaml"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.idx} \
        --readFilesIn {input.r1} {input.r2} \
        --readFilesCommand zcat \
        --outFileNamePrefix results/with_ref_annot/{wildcards.sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
        """
