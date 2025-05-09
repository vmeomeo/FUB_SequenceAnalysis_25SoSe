# Test
# print(samples.at["ERR024604_tiny", "fq2"]) 

rule create_index:
    input: 
        lambda wc: f"{config['ref']}"
    output:
        expand(
            "{output_dir}/ref_idx/reference_index.{ext}",
            output_dir=config["output_dir_path"],
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    threads: config["threads"]
    conda:
        "../env/bowtie2.yaml"
    shell:
        "bowtie2-build -f {input} {config[output_dir_path]}/ref_idx/reference_index"

print(f"Running Bowtie2 with output: {config['output_dir_path']}/sam/{{sample}}.sam")

rule align:
    input: 
        fq1 = lambda wc: samples.at[wc.sample, "fq1"],
        fq2 = lambda wc: samples.at[wc.sample, "fq2"],
        index = expand(
            "{output_dir}/ref_idx/reference_index.{ext}",
            output_dir=config["output_dir_path"],
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        )
    output:
        f"{config['output_dir_path']}/sam/{{sample}}.sam"
    threads: config["threads"]
    conda:
        "../env/bowtie2.yaml"
    shell: 
        "bowtie2 --no-unal -t --threads {config[bowtie_threads]} "
        "--np {config[N_penalty]} "
        "--rfg {config[ref_gap_penalty]} " # Reason caused the errors --> solved
        "-x {config[output_dir_path]}/ref_idx/reference_index "
        "-1 {input.fq1} -2 {input.fq2} -S {output[0]}"
