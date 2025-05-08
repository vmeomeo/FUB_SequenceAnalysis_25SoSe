### How to run snakemake?

1. Go to directory 'workflow', which has Snakefile inside
```cd your_path/workflow```

2. Dry-run snakemake to check errors
```snakemake -s Snakefile -n```

3. If there are no errors, run (you can change the number of cores)
```snakemake -s Snakefile --cores 4```