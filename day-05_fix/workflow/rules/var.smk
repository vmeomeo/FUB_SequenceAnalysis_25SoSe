rule var:
    input:
        alignment = f"{config['output_dir_path']}/alignment/all_samples_aligned.fasta"
    output:
        variability_txt = f"{config['output_dir_path']}/variability/all_samples_windows_var.txt",
        plot_png = f"{config['output_dir_path']}/variability/all_samples_var_plot.png"
    params:
        window_size = config['window_size'],
        step_size = config['step_size']
    conda:
        "../env/var.yaml"
    log:
        f"{config['output_dir_path']}/variability/all_samples_variability.log"
    shell:
        """
        workflow/scripts/calc_variability.py {input.alignment} {params.window_size} {output.variability_txt} {output.plot_png} {params.step_size} > {log} 2>&1
        """
