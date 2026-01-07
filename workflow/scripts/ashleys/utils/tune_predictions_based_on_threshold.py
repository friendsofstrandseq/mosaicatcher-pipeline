import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t")
df.loc[df["probability"] >= float(snakemake.config["ashleys_threshold"]), "prediction"] = 1
df.loc[df["probability"] < float(snakemake.config["ashleys_threshold"]), "prediction"] = 0
df.sort_values(by="cell").to_csv(snakemake.output[0], sep="\t", index=False)