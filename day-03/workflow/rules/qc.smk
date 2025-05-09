rule fastqc:
    input:
        fq1 = lambda wc: samples.at[wc.sample, "fq1"],
        fq2 = lambda wc: samples.at[wc.sample, "fq2"]
    output: 
    params:
    threads: 8
    conda: "../env/fastqc.yaml"
    wrapper: 
        "fastqc {input.fq1} {input.fq2}"

        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN