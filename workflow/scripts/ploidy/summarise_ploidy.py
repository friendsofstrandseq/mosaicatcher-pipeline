import pandas as pd

df = pd.read_csv(snakemake.input.ploidy, sep="\t")
df = df.groupby("#chrom")["ploidy_estimate"].describe()
df.to_csv(snakemake.output.summary, sep="\t")
