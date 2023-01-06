import pandas as pd
import numpy as np

ploidy_detailed = pd.read_csv(snakemake.input[0], sep="\t")
ploidy_detailed = ploidy_detailed.loc[ploidy_detailed["#chrom"] != "genome"]

m_f = "M"

x_mean = ploidy_detailed.loc[ploidy_detailed["#chrom"] == "chrX", "ploidy_estimate"].mean()
if x_mean > 0:
    m_f = "F" if x_mean >= 2 else "M"

ploidy_detailed["sex"] = m_f
ploidy_detailed["start"] = ploidy_detailed["start"] + 1
ploidy_detailed.loc[ploidy_detailed["ploidy_estimate"] > 2, "ploidy_estimate"] = 2
ploidy_detailed = ploidy_detailed[["#chrom", "start", "end", "sex", "ploidy_estimate"]]
ploidy_detailed.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
