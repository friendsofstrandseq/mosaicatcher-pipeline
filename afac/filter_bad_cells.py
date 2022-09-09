import subprocess
import pandas as pd
import sys, os

# snakemake_log = open(snakemake.log[0], "w")

# Prepare header of info files
# subprocess.call("grep '^#' {} > {}".format(snakemake.input.info_raw, snakemake.output.info), shell=True)
# subprocess.call("grep '^#' {} > {}".format(snakemake.input.info_raw, snakemake.output.info_removed), shell=True)
snakemake_input_info_raw = ".tests/output_CHR17_bam_singularity/counts/RPE-BM510/RPE-BM510.info_raw"
snakemake_output_info = ".tests/output_CHR17_bam_singularity/counts/RPE-BM510/RPE-BM510.info"
snakemake_output_info_removed = ".tests/output_CHR17_bam_singularity/counts/RPE-BM510/RPE-BM510.info_rm"


subprocess.call("grep '^#' {} > {}".format(snakemake_input_info_raw, snakemake_output_info), shell=True)
subprocess.call("grep '^#' {} > {}".format(snakemake_input_info_raw, snakemake_output_info_removed), shell=True)

# Read mosaic count info
df = pd.read_csv(snakemake_input_info_raw, skiprows=13, sep="\t")
df["pass1"] = df["pass1"].astype(int)

labels_path = ".tests/output_CHR17_bam_singularity/config/RPE-BM510/labels.tsv"
# labels_path = snakemake.input.labels
labels = pd.read_csv(labels_path, sep="\t")


# snakemake_log.write(labels.to_str())

# b_ashleys = "ENABLED" if snakemake.config["ashleys_pipeline"] is True else "DISABLED"
# b_old = "ENABLED" if snakemake.config["input_old_behavior"] is True else "DISABLED"

b_old = "ENABLED"

# snakemake_log.write("ASHLEYS preprocessing module: {}".format(b_ashleys))
# snakemake_log.write("input_old_behavior parameter: {}".format(b_old))
# snakemake_log.write("Computing intersection between lists ...")

# IF BOTH MOSAIC INFO FILE & LABELS DF ARE AVAILABLE + SAME SIZE
# if labels.shape[0] == df.shape[0]:

cells_to_keep_labels = labels.loc[labels["prediction"] == 1]["cell"].str.replace(".sort.mdup.bam", "").sort_values().tolist()
cells_to_keep_mosaic = df.loc[df["pass1"] == 1]["cell"].unique().tolist()
cells_to_keep = list(sorted(list(set(cells_to_keep_labels).intersection(cells_to_keep_mosaic))))

# else:
# CATCH ERROR IF DIFFERENT SIZES AND CONFIG ENABLED
# if (snakemake.config["ashleys_pipeline"] is True) or (snakemake.config["input_old_behavior"] is True):
#     os.exit("Dataframes do not have the same dimensions:")
#     os.exit("mosaic info: {} ; labels: {}".format(str(df.shape[0]), str(labels.shape[0])))

# # ELSE NORMAL MODE
# else:
# snakemake_log.write("Standard mode using only 'mosaic count info' file")
# cells_to_keep = df.loc[df["pass1"] == 1]["cell"].unique().tolist()


# cells_to_keep = labels.loc[labels["prediction"] == 1]["cell"].str.replace(".sort.mdup.bam", "").tolist()
df_kept = df.loc[df["cell"].isin(cells_to_keep)]
df_removed = df.loc[~df["cell"].isin(cells_to_keep)]

print(cells_to_keep)

# snakemake_log.write("List of cells kept: ")
# for cell in sorted(cells_to_keep):
#     snakemake_log.write("- {cell}".format(cell=cell))

# snakemake_log.write("List of cells removed:")
# for cell in sorted(df_removed["cell"].values.tolist()):
#     snakemake_log.write("- {cell}".format(cell=cell))


# df_counts = pd.read_csv(snakemake.input.counts_sort, compression="gzip", sep="\t")
# df_counts = df_counts.loc[df_counts["cell"].isin(cells_to_keep)]
# # print(df_counts)/
# df_counts.to_csv(snakemake.output.counts, compression="gzip", sep="\t", index=False)

# df_kept.to_csv(snakemake.output.info, index=False, sep="\t", mode="a")
# df_removed.to_csv(snakemake.output.info_removed, index=False, sep="\t", mode="a")
