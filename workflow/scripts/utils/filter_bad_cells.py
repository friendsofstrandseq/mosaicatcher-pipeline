import subprocess
import pandas as pd

subprocess.call("grep '^#' {} > {}".format(snakemake.input.info_raw, snakemake.output.info), shell=True)
subprocess.call("grep '^#' {} > {}".format(snakemake.input.info_raw, snakemake.output.info_removed), shell=True)

df = pd.read_csv(snakemake.input.info_raw, skiprows=13, sep="\t")
df["pass1"] = df["pass1"].astype(int)
if snakemake.config["ashleys_pipeline"] is False:
    print("Ashleys pipeline NOT RUNNED, filtering cells that cannot be used for segmentation:")

    df_kept = df.loc[df["pass1"] == 1]
    df_removed = df.loc[df["pass1"] == 0]
    print(df_removed["Cell"].unique().tolist())
    cells_to_keep = df_removed["Cell"].unique().tolist()
elif snakemake.config["ashleys_pipeline"] is True:
    # if os.path.isfile(labels_path) is True:
    print("Ashleys pipeline RUNNED, filtering cells that cannot be used for segmentation based on ashleys-qc predictions:")
    labels_path = "{folder}/{sample}/cell_selection/labels.tsv".format(
        folder=snakemake.config["input_bam_location"], sample=snakemake.wildcards.sample
    )
    # labels_path = snakemake.input.labels
    labels = pd.read_csv(labels_path, sep="\t")
    print(labels)
    cells_to_keep = labels.loc[labels["prediction"] == 1]["cell"].str.replace(".sort.mdup.bam", "").tolist()
    df_kept = df.loc[df["cell"].isin(cells_to_keep)]
    df_removed = df.loc[~df["cell"].isin(cells_to_keep)]


print(sorted(cells_to_keep))

df_counts = pd.read_csv(snakemake.input.counts_sort, compression="gzip", sep="\t")
df_counts = df_counts.loc[df_counts["cell"].isin(cells_to_keep)]
df_counts.to_csv(snakemake.output.counts, compression="gzip", sep="\t", index=False)

df_kept.to_csv(snakemake.output.info, index=False, sep="\t", mode="a")
df_removed.to_csv(snakemake.output.info_removed, index=False, sep="\t", mode="a")
