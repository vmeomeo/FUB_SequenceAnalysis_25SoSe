rule decompress_fasta:
    input:
        "resources/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
    output:
        "resources/Homo_sapiens.GRCh38.dna.chromosome.22.fa"
    shell:
        "gunzip -c {input} > {output}"

rule decompress_gtf:
    input:
        "resources/Homo_sapiens.GRCh38.114.chr.gtf.gz"
    output:
        "resources/Homo_sapiens.GRCh38.114.chr.gtf"
    shell:
        "gunzip -c {input} > {output}"
rule genome_indexing:
    input:
        fa = "resources/Homo_sapiens.GRCh38.dna.chromosome.22.fa",
        gtf = "resources/Homo_sapiens.GRCh38.114.chr.gtf"
    output:
        directory("resources/index")
    threads: 4
    conda:
        "../env/star.yaml"
    log:
        "results/indexing/logs/STAR_genome_indexing.log"
    shell:
        """
        mkdir -p {output}
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
        idx = directory("resources/index"),
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
