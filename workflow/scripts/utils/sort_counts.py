import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t", compression="gzip")
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)
chroms = ["chr{}".format(str(c)) for c in list(range(1, 23))] + ["chrX", "chrY"]
df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
df.sort_values(by=["cell", "chrom", "start", "end"]).to_csv(snakemake.output[0], compression="gzip", sep="\t", index=False)
