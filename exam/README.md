### How to run snakemake?

1. Go to directory outside 'workflow' 1 level

2. Dry-run snakemake to check errors
```snakemake --use-conda --cores 16 -n```

3. If there are no errors, run (you can change the number of cores)
```snakemake --use-conda --cores 16```