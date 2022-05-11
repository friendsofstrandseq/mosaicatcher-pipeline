import pandas as pd
import sys

df_counts = pd.read_csv(sys.argv[1], compression="gzip", sep="\t")
sample_cells = df_counts["cell"].unique().tolist()[:2]
df_counts_lite = df_counts.loc[df_counts["cell"].isin(sample_cells)]
df_counts_lite.to_csv(sys.argv[3], sep="\t", compression="gzip", index=False)

info_header = "".join([line for line in open(sys.argv[2], "r").readlines() if line[0] == "#"])
with open(sys.argv[4], "w") as w:
    w.write(info_header)

df_info = pd.read_csv(sys.argv[2], sep="\t", skiprows=13)
df_info_lite = df_info.loc[df_info["cell"].isin(sample_cells)]
df_info_lite.to_csv(sys.argv[4], sep="\t", mode="a", index=False)
