import os

# Define the parent directory containing sample subdirectories
parent_dir = "results/assembly"

# Loop through each subdirectory
for sample in os.listdir(parent_dir):
    sample_dir = os.path.join(parent_dir, sample)
    if os.path.isdir(sample_dir):
        for fname in os.listdir(sample_dir):
            if fname.startswith("assembly."):
                old_path = os.path.join(sample_dir, fname)
                new_fname = fname.replace("assembly", sample, 1)
                new_path = os.path.join(sample_dir, new_fname)
                print(f"Renaming: {old_path} -> {new_path}")
                os.rename(old_path, new_path)