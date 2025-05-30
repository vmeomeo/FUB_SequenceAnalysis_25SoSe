# scripts/aggregate_idxstats.py

import pandas as pd
from pathlib import Path

dfs = []
for path in snakemake.input.idxstats:
    sample = Path(path).stem
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["ref", "length", "mapped", "unmapped"])
    df = df[["ref", "length", "mapped"]].copy()
    df.rename(columns={"mapped": sample}, inplace=True)
    dfs.append(df.set_index(["ref", "length"]))

merged = pd.concat(dfs, axis=1).reset_index()
merged.to_csv(snakemake.output[0], sep="\t", index=False)
