import pandas as pd
import os
from snakemake import input, output, log

log_file = open(log[0], "w")

def concat_tables(files, source):
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file, sep="\t")
            df.insert(0, "Sample", os.path.basename(file).split(".")[0])
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}", file=log_file)
    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        return combined
    else:
        return pd.DataFrame()

# Load all modules
mlst_df = concat_tables(snakemake.input.mlst, "mlst")
card_df = concat_tables(snakemake.input.card, "card")
vfdb_df = concat_tables(snakemake.input.vfdb, "vfdb")
mob_df  = concat_tables(snakemake.input.mob, "mob")

# Write to Excel with multiple sheets
with pd.ExcelWriter(snakemake.output.summary) as writer:
    mlst_df.to_excel(writer, sheet_name="MLST", index=False)
    card_df.to_excel(writer, sheet_name="Resistance_CARD", index=False)
    vfdb_df.to_excel(writer, sheet_name="Virulence_VFDB", index=False)
    mob_df.to_excel(writer, sheet_name="Plasmids", index=False)
