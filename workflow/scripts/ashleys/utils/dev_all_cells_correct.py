import pandas as pd
import os
import glob

# Handle Namedlist from Snakemake v9 - select_ashleys_labels returns expand() which is always a list/Namedlist
# Even with one element, we need to index it: snakemake.input.folder[0]
try:
    input_file = str(snakemake.input.folder[0])
except (TypeError, IndexError, AttributeError):
    input_file = str(snakemake.input.folder)

df = pd.read_csv(input_file, sep="\t")
df["prediction"] = 1
df["probability"] = 1
df.loc[df["cell"].str.contains("BM510x04_PE20305"), "prediction"] = 0
df.loc[df["cell"].str.contains("BM510x04_PE20312"), "prediction"] = 0

# Ensure all cells from FASTQ directory are included (important for downsampled datasets)
fastq_dir = f"{snakemake.config['data_location']}/{snakemake.wildcards.sample}/fastq"
if os.path.isdir(fastq_dir):
    # Get unique cell names from fastq files (remove .1.fastq.gz and .2.fastq.gz suffixes)
    fastq_files = sorted(glob.glob(f"{fastq_dir}/*.fastq.gz"))
    all_cells = sorted(set([os.path.basename(f).rsplit('.', 3)[0] + ".sort.mdup.bam" for f in fastq_files]))

    # Add missing cells with prediction=1, probability=1 (default for light data)
    missing_cells = [cell for cell in all_cells if cell not in df["cell"].values]
    if missing_cells:
        missing_df = pd.DataFrame({
            "cell": missing_cells,
            "prediction": [1] * len(missing_cells),
            "probability": [1] * len(missing_cells)
        })
        df = pd.concat([df, missing_df], ignore_index=True)

# if "BM510" in snakemake.wildcards.sample and snakemake.config["use_light_data"] is True:
#     df = pd.concat([df, pd.DataFrame([{"cell": "BM510x04_PE20320.sort.mdup.bam", "prediction": 0, "probability": 0}])])
df = df.sort_values(by="cell")
df.to_csv(snakemake.output.folder, sep="\t", index=False)
