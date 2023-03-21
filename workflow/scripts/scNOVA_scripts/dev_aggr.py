import pandas as pd
import os

l = list(snakemake.input)

l_df = list()

l_df.append(
    pd.read_csv(
        sorted(list(l))[0],
        sep="\t",
        names=[
            "'chr'",
            "'start'",
            "'end'",
            "Feature",
            "'{file}'".format(file=os.path.basename(sorted(list(l))[0])).replace(".tab", ".bam"),
        ],
    )[["'chr'", "'start'", "'end'"]]
)

l_df.extend(
    [
        pd.read_csv(
            file,
            sep="\t",
            names=[
                "'chr'",
                "'start'",
                "'end'",
                "Feature",
                "'{file}'".format(file=os.path.basename(file)).replace(".tab", ".bam"),
            ],
        )[["'{file}'".format(file=os.path.basename(file)).replace(".tab", ".bam")]]
        for file in sorted(list(l))
    ]
)
# print(len(l_df))

pd.concat(l_df, axis=1).to_csv(snakemake.output.tab, sep="\t", index=False)
