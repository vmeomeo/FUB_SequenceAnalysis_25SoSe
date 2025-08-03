import pandas as pd
import os

log_file = open(snakemake.log[0], "w")

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

with pd.ExcelWriter(snakemake.output.summary) as writer:
    if "mlst" in snakemake.input:
        mlst_df = concat_tables(snakemake.input.mlst, "mlst")
        mlst_df.to_excel(writer, sheet_name="MLST", index=False)

    if "card" in snakemake.input:
        card_df = concat_tables(snakemake.input.card, "card")
        card_df.to_excel(writer, sheet_name="Resistance_CARD", index=False)

    if "vfdb" in snakemake.input:
        vfdb_df = concat_tables(snakemake.input.vfdb, "vfdb")
        vfdb_df.to_excel(writer, sheet_name="Virulence_VFDB", index=False)

    if "mob" in snakemake.input:
        mob_df = concat_tables(snakemake.input.mob, "mob")
        mob_df.to_excel(writer, sheet_name="Plasmids", index=False)
