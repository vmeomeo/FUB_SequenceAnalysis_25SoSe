rule quast:
    input: 
        assembly = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta"
    output:
        assembly_qc = directory(f"{config['output_dir_path']}/qc/quast/{{sample}}")
    conda:
        "../env/quast.yaml"
    threads: workflow.cores
    log:
        f"{config['output_dir_path']}/qc/quast/logs/quast_{{sample}}.log"
    shell:
        """
        quast.py {input.assembly} \
                 --gene-finding \
                 -o {output.assembly_qc} \
                 -t {threads} \
                 > {log} 2>&1
        """

rule download_busco_db:
    output:
        touch(f"{config['busco_db_path']}/.downloaded")
    log:
        f"{config['output_dir_path']}/logs/busco/download_busco_db.log"
    conda:
        "../env/busco.yaml"
    params:
        lineage = config["bacteria_odb10"]
    shell:
        """
        busco --download {params.lineage} --download_path {config[busco_db_path]} &> {log}
        touch {output}
        """


rule busco:
    input:
        fasta = f"{config['output_dir_path']}/assembly/{{sample}}/{{sample}}.fasta",
        db_marker = f"{config['busco_db_path']}/.downloaded"
    output:
        summary = f"{config['output_dir_path']}/qc/busco/{{sample}}/{{sample}}/short_summary.specific.bacteria_odb10.{{sample}}.txt",
        busco_dirs = directory(f"{config['output_dir_path']}/qc/busco/{{sample}}/{{sample}}")
    conda:
        "../env/busco.yaml"
    params:
        lineage = config["bacteria_odb10"],
        outdir = config['output_dir_path']
    threads:
        workflow.cores
    shell:
        """
        busco -i {input.fasta} \
              -o {wildcards.sample} \
              -l {params.lineage} \
              -m genome \
              -f \
              --out_path {params.outdir}/qc/busco/{wildcards.sample} \
              --cpu {threads}
        """


rule assembly_multiqc:
    input:
        quast_dirs = expand(f"{config['output_dir_path']}/qc/quast/{{sample}}", sample=get_included_samples()),
        # busco_summaries = expand(
        #     f"{config['output_dir_path']}/qc/busco/{{sample}}/{{sample}}/short_summary.specific.bacteria_odb10.{{sample}}.txt",
        #     sample=samples.index
        # )
        busco_dirs = expand(
            f"{config['output_dir_path']}/qc/busco/{{sample}}/{{sample}}", sample=get_included_samples()
        )
    output:
        html = f"{config['output_dir_path']}/qc/multiqc_assembly_report.html"
    conda:
        "../env/multiqc.yaml"
    threads: 1
    params:
        outdir = f"{config['output_dir_path']}/qc"
    log:
        f"{config['output_dir_path']}/qc/logs/multiqc_assembly.log"
    shell:
        """
        multiqc {input.quast_dirs} {input.busco_dirs} -o {params.outdir} -n multiqc_assembly_report > {log} 2>&1
        """
        
        # rm -f {params.outdir}/multiqc_assembly_report*.html
        # multiqc {input.quast_dirs} {input.busco_summaries} -o {params.outdir} -n multiqc_assembly_report > {log} 2>&1
        # """
