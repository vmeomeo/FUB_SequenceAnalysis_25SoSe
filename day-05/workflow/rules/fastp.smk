rule fastp:
    input:
        fq = lambda wc: samples.at[wc.sample, f"fq{wc.read_number}"]
    output:
        files = f"{config['output_dir_path']}/cleaned/{{sample}}_{{read_number}}.fastq",
        html = f"{config['output_dir_path']}/qc/{{sample}}_{{read_number}}_fastp.html",
        json = f"{config['output_dir_path']}/qc/{{sample}}_{{read_number}}_fastp.json"
    params:
        config_dir = config['output_dir_path'],
        fastp_flags = config['fastp_flags']
    conda:
        "../env/fastp.yaml"
    log:
        f"{config['output_dir_path']}/qc/logs/fastp_{{sample}}_{{read_number}}.log"
    shell:
        """
        fastp {params.fastp_flags} -i {input.fq} -o {output.files} -j {output.json} -h {output.html} > {log} 2>&1
        """
