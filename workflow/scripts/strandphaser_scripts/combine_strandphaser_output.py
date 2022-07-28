import pandas as pd

l_df = list()
print(snakemake.input.files)
for j, file in enumerate(snakemake.input.files):
    print(j, file)
    tmp_df = pd.read_csv(file, sep="\t")
    print(tmp_df)
    l_df.append(tmp_df)
df = pd.concat(l_df)
df.to_csv(snakemake.output[0], sep="\t", index=False)
