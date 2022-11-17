import pandas as pd
import os
df = pd.read_csv(snakemake.input[0], sep="\t")
df["callset"] = df["callset"].apply(lambda r: os.path.basename(r).replace(".tsv", ""))
df = df.set_index("callset")
df = df.fillna(0).T.reset_index()
pd.options.display.float_format = '{:,.1f}'.format
df_out = write_to_html_file(df, snakemake.wildcards.sample)
with open(snakemake.output.html, 'w') as o:
    o.write(df_out)