# rule unicycler_short:
#     # input:
#     #     r1 = f"{config['output_dir_path']}/cleaned/{{sample}}_r1.fastq.gz",
#     #     r2 = f"{config['output_dir_path']}/cleaned/{{sample}}_r2.fastq.gz"
#     input:
#         r1 = lambda wc: (
#             f"{config['output_dir_path']}/decontaminated/{wc.sample}_r1.fastq.gz"
#             if config.get("contaminant_fasta", "")
#             else f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz"
#         ),
#         r2 = lambda wc: (
#             f"{config['output_dir_path']}/decontaminated/{wc.sample}_r2.fastq.gz"
#             if config.get("contaminant_fasta", "")
#             else f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz"
#         )
#     output:
#         f"{config['output_dir_path']}/assembly/{{sample}}/assembly.fasta"
#     log:
#         f"{config['output_dir_path']}/logs/unicycler/{{sample}}.log"
#     threads: workflow.cores
#     conda:
#         "../env/unicycler.yaml"
#     params:
#         outdir = config['output_dir_path']
#     shell:
#         "unicycler -1 {input.r1} -2 {input.r2} -o {params.outdir}/assembly/{wildcards.sample} "
#         "-t {threads} &> {log}"

# rule unicycler:
#     input:
#         r1 = lambda wc: (
#             f"{config['output_dir_path']}/decontaminated/{wc.sample}_r1.fastq.gz"
#             if config.get("contaminant_fasta", "")
#             else f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz"
#         ),
#         r2 = lambda wc: (
#             f"{config['output_dir_path']}/decontaminated/{wc.sample}_r2.fastq.gz"
#             if config.get("contaminant_fasta", "")
#             else f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz"
#         ),
#         l = lambda wc: (f"{config['output_dir_path']}/filtered/{wc.sample}_filtered.fastq")
#     output:
#         f"{config['output_dir_path']}/assembly/{{sample}}/assembly.fasta"
#     log:
#         f"{config['output_dir_path']}/logs/unicycler/{{sample}}.log"
#     threads: workflow.cores
#     conda:
#         "../env/unicycler.yaml"
#     params:
#         outdir = config['output_dir_path']
#     shell:
#         "unicycler -1 {input.r1} -2 {input.r2} -l {input.l} -o {params.outdir}/assembly/{wildcards.sample} "
#         "-t {threads} &> {log}"

# rule change_assembly_name_unicycler:
#     input:
#         assembly = f"{config['output_dir_path']}/assembly/{{sample}}/assembly.fasta"
#     output:
#         assembly = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta"
#     script:
#         "../scripts/change_assembly_name_fasta.py"

rule unicycler:
    input:
        r1 = lambda wc: (
            f"{config['output_dir_path']}/decontaminated/{wc.sample}_r1.fastq.gz"
            if config.get("contaminant_fasta", "")
            else f"{config['output_dir_path']}/cleaned/{wc.sample}_r1.fastq.gz"
        ),
        r2 = lambda wc: (
            f"{config['output_dir_path']}/decontaminated/{wc.sample}_r2.fastq.gz"
            if config.get("contaminant_fasta", "")
            else f"{config['output_dir_path']}/cleaned/{wc.sample}_r2.fastq.gz"
        ),
        l = lambda wc: f"{config['output_dir_path']}/filtered/{wc.sample}_filtered.fastq"
    output:
        renamed = directory(f"{config['output_dir_path']}/assembly/{{sample}}")
    log:
        f"{config['output_dir_path']}/logs/unicycler/{{sample}}.log"
    threads: workflow.cores
    conda:
        "../env/unicycler.yaml"
    params:
        outdir = f"{config['output_dir_path']}/assembly"
    run:
        import os, glob, shutil

        sample = wildcards.sample
        tmpdir = f"{params.outdir}/{sample}/_unicycler_tmp"
        finaldir = f"{params.outdir}/{sample}"

        shell(
            f"unicycler -1 {input.r1} -2 {input.r2} -l {input.l} "
            f"-o {tmpdir} -t {threads} &> {log}"
        )

        # Ensure final output dir exists
        os.makedirs(finaldir, exist_ok=True)

        # Rename all output files
        for f in glob.glob(f"{tmpdir}/*"):
            basename = os.path.basename(f)
            if basename.startswith("assembly"):
                ext = basename.replace("assembly", "")
                new_name = f"{sample}{ext}"
                shutil.move(f, os.path.join(finaldir, new_name))
            else:
                shutil.move(f, os.path.join(finaldir, basename))