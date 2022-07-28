import pandas as pd

l_df = list()
for j, file in enumerate(snakemake.input[0]):
    tmp_df = pd.read_csv(file, sep="\t")
    print(tmp_df)
    l_df.append(tmp_df)
df = pd.concat(l_df)
df.to_csv(snakemake.output[0], sep="\t", index=False)
