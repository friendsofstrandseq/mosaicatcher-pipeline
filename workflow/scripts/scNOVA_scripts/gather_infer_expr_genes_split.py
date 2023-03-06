import pandas as pd

pd.concat([pd.read_csv(e) for e in sorted(list(snakemake.input))])[["prob1", "prob2"]].reset_index().to_csv(
    snakemake.output[0], index=False
)
