import pandas as pd
from pathlib import Path
import sys

def process_inputs(input_files, output_xlsx):
    writer = pd.ExcelWriter(output_xlsx, engine='openpyxl')
    
    # Group files by type
    file_groups = {}
    for f in input_files:
        path = Path(f)
        if 'mlst' in str(path):
            file_groups.setdefault('mlst', []).append(f)
        elif 'card' in str(path):
            file_groups.setdefault('card', []).append(f)
        elif 'vfdb' in str(path):
            file_groups.setdefault('vfdb', []).append(f)
        elif 'plasmids' in str(path):
            file_groups.setdefault('mob', []).append(f)
    
    # Process each group
    for sheet_name, files in file_groups.items():
        dfs = []
        for f in files:
            try:
                df = pd.read_csv(f, sep='\t')
                sample = Path(f).stem.split('.')[0]
                df['Sample'] = sample  # Add sample identifier
                dfs.append(df)
            except Exception as e:
                print(f"Error processing {f}: {e}", file=sys.stderr)
        
        if dfs:
            pd.concat(dfs).to_excel(writer, sheet_name=sheet_name[:31], index=False)
    
    writer.close()

if __name__ == '__main__':
    input_files = snakemake.input  # List of all input files
    output_xlsx = snakemake.output.summary
    Path(output_xlsx).parent.mkdir(parents=True, exist_ok=True)
    process_inputs(input_files, output_xlsx)