import os, sys
import pandas as pd
import scipy

# LOAD MOSAIC COUNTS INFO
info = snakemake.input.info
info_df = pd.read_csv(info, sep="\t", skiprows=13)
info_df["cell"] = info_df["cell"] + ".sort.mdup.bam"

print(f"Info path: {info}")

# Load ashleys predictions
labels_path = snakemake.input.labels[0]

print(f"Labels path: {labels_path}")
print(f"Labels path type: {type(labels_path)}")
labels = pd.read_csv(labels_path, sep="\t").sort_values(by="cell")
labels["sample"] = snakemake.wildcards.sample
print(labels)
# Neg control with no reads
info_df_neg_control = info_df.loc[info_df["mapped"] == 0]
print(info_df_neg_control)
print(f"Shape of info_df_neg_control: {info_df_neg_control.shape}")
if info_df_neg_control.shape[0] > 0:
    info_df_neg_control["prediction"] = 0
    info_df_neg_control["probability"] = 0
    info_df_neg_control["sample"] = snakemake.wildcards.sample
    labels_neg_control = info_df_neg_control[
        ["cell", "prediction", "probability", "sample"]
    ]
    print(labels_neg_control)
    labels = pd.concat([labels, labels_neg_control]).sort_values(by="cell")
print(labels)

# Retrieve correct prediction
labels_corrected = labels.loc[labels["prediction"] == 1].sort_values(by="cell")

# Merge counts & labels
labels_corrected = pd.merge(info_df[["cell", "good"]], labels_corrected, on="cell")
# Compute z-score on reads nb
labels_corrected["z_score"] = scipy.stats.zscore(labels_corrected["good"])

print(labels_corrected)

# Output, Correct outliers predictions & proba
z_score_cutoff = 5
labels_corrected.loc[labels_corrected["z_score"] >= z_score_cutoff, "cell"].to_csv(
    snakemake.output.bypass_cell, index=False, sep="\t"
)
print(labels_corrected)
if labels_corrected.loc[labels_corrected["z_score"] >= z_score_cutoff].shape[0] > 0:
    labels_corrected.loc[
        labels_corrected["z_score"] >= z_score_cutoff, "new_prediction"
    ] = 0
    labels_corrected.loc[
        labels_corrected["z_score"] >= z_score_cutoff, "new_probability"
    ] = 0
    labels_corrected.loc[
        labels_corrected["z_score"] < z_score_cutoff, "new_probability"
    ] = labels_corrected.loc[
        labels_corrected["z_score"] < z_score_cutoff, "probability"
    ]
    labels_corrected["new_prediction"] = labels_corrected["new_prediction"].fillna(1)
    labels_corrected["new_prediction"] = labels_corrected["new_prediction"].astype(int)

    # Back to full dataframe
    labels.loc[
        labels["cell"].isin(labels_corrected.cell.values.tolist()), "prediction"
    ] = labels_corrected.new_prediction.values.tolist()
    labels.loc[
        labels["cell"].isin(labels_corrected.cell.values.tolist()), "probability"
    ] = labels_corrected.new_probability.values.tolist()
else:
    print("All cells were discarded")
    # labels = labels_corrected

print(labels)

# Output
labels.sort_values(by="cell").to_csv(
    snakemake.output.labels_corrected, index=False, sep="\t"
)
