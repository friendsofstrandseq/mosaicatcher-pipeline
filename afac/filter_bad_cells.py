import subprocess
import pandas as pd


config = {"ashleys_pipeline": True}

info_raw = ".tests/output_T2T/counts/RPE-BM510/RPE-BM510.info_raw"
labels_path = ".tests/data_T2T/RPE-BM510/cell_selection/labels.tsv"


df = pd.read_csv(info_raw, skiprows=13, sep="\t")
df["pass1"] = df["pass1"].astype(int)
if config["ashleys_pipeline"] is False:
    print("Ashleys pipeline NOT RUNNED, filtering cells that cannot be used for segmentation:")

    df_kept = df.loc[df["pass1"] == 1]
    df_removed = df.loc[df["pass1"] == 0]
    # print(df_removed["Cell"].unique().tolist())
    cells_to_keep = df_kept["cell"].unique().tolist()
elif config["ashleys_pipeline"] is True:
    # if os.path.isfile(labels_path) is True:
    print("Ashleys pipeline RUNNED, filtering cells that cannot be used for segmentation based on ashleys-qc predictions:")
    # labels_path = snakemake.input.labels
    labels = pd.read_csv(labels_path, sep="\t")

    cells_to_keep_ashleys = labels.loc[labels["prediction"] == 1]["cell"].str.replace(".sort.mdup.bam", "").sort_values().tolist()
    cells_to_keep_mosaic = df.loc[df["pass1"] == 1]["cell"].unique().tolist()
    cells_to_keep = list(sorted(list(set(cells_to_keep_ashleys).intersection(cells_to_keep_mosaic))))

    df_kept = df.loc[df["cell"].isin(cells_to_keep)]
    df_removed = df.loc[~df["cell"].isin(cells_to_keep)]


# print(sorted(cells_to_keep))

# df_counts = pd.read_csv(snakemake.input.counts_sort, compression="gzip", sep="\t")
# df_counts = df_counts.loc[df_counts["cell"].isin(cells_to_keep)]
# df_counts.to_csv(snakemake.output.counts, compression="gzip", sep="\t", index=False)

# df_kept.to_csv(snakemake.output.info, index=False, sep="\t", mode="a")
# df_removed.to_csv(snakemake.output.info_removed, index=False, sep="\t", mode="a")
