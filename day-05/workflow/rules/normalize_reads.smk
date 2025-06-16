rule normalize_reads:
    input:
        r1 = f"{config['output_dir_path']}/decontaminated/{{sample}}_1.fastq",
        r2 = f"{config['output_dir_path']}/decontaminated/{{sample}}_2.fastq"
    output:
        r1_norm = f"{config['output_dir_path']}/normalized/{{sample}}_1.fastq",
        r2_norm = f"{config['output_dir_path']}/normalized/{{sample}}_2.fastq"
    conda:
        "../env/bbnorm.yaml"
    threads: 4
    log:
        f"{config['output_dir_path']}/normalized/logs/{{sample}}.log"
    params:
        flags = config.get("bbnorm_flags", "target=40 min=5")  # customize if needed
    shell:
        """
        bbnorm.sh in1={input.r1} in2={input.r2} \
                  out1={output.r1_norm} out2={output.r2_norm} \
                  threads={threads} {params.flags} > {log} 2>&1
        """
