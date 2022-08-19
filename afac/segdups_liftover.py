import pandas as pd

df_segdups_hg38_full = pd.read_csv(
    "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/data/segdups/segDups_hg38_UCSCtrack.bed.gz", compression="gzip", sep="\t"
)
col_list = list(df_segdups_hg38_full.columns)
df_segdups_hg38_full = df_segdups_hg38_full.rename(
    {"chrom": "chrom_hg38", "chromStart": "chromStart_hg38", "chromEnd": "chromEnd_hg38"}, axis=1
)

df_segdups_new_lite = pd.read_csv(
    "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/data/segdups/segDups_T2T_lite.bed",
    sep="\t",
    names=["chromStart", "chromEnd", "name"],
)
df_segdups_new_lite.index.name = "chrom"
df_segdups_new_lite = df_segdups_new_lite.reset_index()

merge_df = pd.merge(df_segdups_hg38_full, df_segdups_new_lite, on="name")
merge_df = merge_df[col_list + ["chrom_hg38", "chromStart_hg38", "chromEnd_hg38"]]

print(merge_df)
merge_df.to_csv(
    "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/data/segdups/segDups_T2T_UCSCtrack.bed.gz",
    index=False,
    sep="\t",
    compression="gzip",
)
