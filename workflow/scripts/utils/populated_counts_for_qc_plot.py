import pandas as pd

# Read 200kb bins file
binbed = pd.read_csv(
    snakemake.input.bin_bed,
    # "../../../../mosaicatcher-update/workflow/data/bin_200kb_all.bed",
    sep="\t",
    names=["chrom", "start", "end", "bin_id"],
)
binbed["ID"] = (
    binbed["chrom"]
    + "_"
    + binbed["start"].astype(str)
    + "_"
    + binbed["end"].astype(str)
)

# Turn chrom into categorical
binbed["chrom"] = pd.Categorical(
    binbed["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)

# Sort & filter out chrY #TMP / can be changed
binbed = binbed.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
binbed["w"], binbed["c"], binbed["class"] = 0, 0, None


# Read SV file
# df = pd.read_csv("../../../../mosaicatcher-update/.tests/data_CHR17/RPE-BM510/counts/RPE-BM510.txt.raw.gz", sep="\t")
df = pd.read_csv(snakemake.input.counts, sep="\t")
df["ID"] = df["chrom"] + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)

## Populate counts df for each cell in order to have all bins represented
l = list()

# Loop over cells
for cell in df.cell.unique().tolist():

    # Outer join to retrieve both real count values from specified chromosome and empty bins
    tmp_df = pd.concat(
        [
            binbed.loc[
                ~binbed["ID"].isin(df.loc[df["cell"] == cell].ID.values.tolist())
            ],
            df.loc[df["cell"] == cell],
        ]
    )

    # Filla cell & sample columns
    tmp_df["cell"] = cell
    tmp_df["sample"] = df.loc[df["cell"] == cell, "sample"].values.tolist()[0]
    l.append(tmp_df)

# Concat list of DF and output
populated_df = pd.concat(l).sort_values(by=["cell", "chrom", "start"])
# populated_df.to_csv("test.txt.gz", compression="gzip", sep="\t", index=False)
populated_df.to_csv(
    snakemake.output.populated_counts, compression="gzip", sep="\t", index=False
)
