import pandas as pd
import yaml

labels = snakemake.input.labels
info_raw = snakemake.input.info_raw
ploidy_summary = snakemake.input.ploidy_summary
single_paired_end_detect = snakemake.input.single_paired_end_detect

single_paired_end_detect_content = open(single_paired_end_detect, "r").readlines()[0]

df_labels = pd.read_csv(labels, sep="\t")
df_labels = df_labels[["cell", "prediction"]]
df_labels = df_labels.rename({"prediction": "hand_labels"}, axis=1)
df_labels["cell"] = df_labels["cell"].str.replace(".sort.mdup.bam", "")

df_info = pd.read_csv(info_raw, skiprows=13, sep="\t")[["cell", "pass1"]]
df_info = df_info.rename({"pass1": "mosaic_cov_pass"}, axis=1)

if df_labels.shape[0] > 0:
    final_df = pd.merge(df_labels, df_info, on="cell")
    final_df.loc[(final_df["hand_labels"] == 1) & (final_df["mosaic_cov_pass"] == 1), "Final_keep"] = 1
    final_df["Final_keep"] = final_df["Final_keep"].fillna(0)

else:
    final_df = df_info
    final_df["Final_keep"] = final_df["mosaic_cov_pass"]
final_df = final_df.rename({"hand_labels": "Ashleys/hand labels"}, axis=1).sort_values(by="cell", ascending=True)


df_ploidy = pd.read_csv(ploidy_summary, sep="\t")[["#chrom", "50%"]]
df_ploidy = df_ploidy.loc[df_ploidy["#chrom"] != "genome"]
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX"]
df_ploidy["#chrom"] = pd.Categorical(df_ploidy["#chrom"], categories=chroms, ordered=True)
df_ploidy = df_ploidy.sort_values(by=["#chrom"]).rename({"#chrom": "chrom", "50%": "ploidy_estimation"}, axis=1)
df_ploidy.loc[df_ploidy["ploidy_estimation"] == 1, "StrandPhaseR_processed"] = 0
df_ploidy["StrandPhaseR_processed"] = df_ploidy["StrandPhaseR_processed"].fillna(1)

with open(snakemake.output.summary, "w") as o:
    o.write("\n==============Library quality summary==============\n")
    o.write("\n")
    o.write(final_df.to_markdown(tablefmt="github", index=False))
    o.write("\n")
    o.write("\n==============Ploidy summary==============\n")
    o.write("\n")
    o.write(df_ploidy.to_markdown(tablefmt="github", index=False))
    o.write("\n")
    o.write("\n==============YAML configuration used==============\n")
    o.write("\n")
    o.write(yaml.dump(snakemake.config))
