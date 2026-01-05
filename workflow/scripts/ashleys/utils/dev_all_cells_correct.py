import pandas as pd

df = pd.read_csv(snakemake.input.folder, sep="\t")
df["prediction"] = 1
df["probability"] = 1
df.loc[df["cell"].str.contains("BM510x04_PE20305"), "prediction"] = 0
df.loc[df["cell"].str.contains("BM510x04_PE20312"), "prediction"] = 0
# if "BM510" in snakemake.wildcards.sample and snakemake.config["use_light_data"] is True:
#     df = pd.concat([df, pd.DataFrame([{"cell": "BM510x04_PE20320.sort.mdup.bam", "prediction": 0, "probability": 0}])])
df = df.sort_values(by="cell")
df.to_csv(snakemake.output.folder, sep="\t", index=False)
