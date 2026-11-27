"""Library size normalisation with strandtools (replaces library_size_normalisation.R).

Median-of-ratios scaling per cell, computed by
strandtools.normalization.library_normalization.
"""

import sys
import types

from tqdm import std as tqdm_std

# strandtools imports tqdm.notebook, which needs ipywidgets: use the console bar.
sys.modules["tqdm.notebook"] = types.SimpleNamespace(tqdm=tqdm_std.tqdm)

import matplotlib

matplotlib.use("Agg")

import pandas as pd
from strandtools.normalization import library_normalization

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

counts = pd.read_csv(snakemake.input.counts, sep="\t", low_memory=False)
counts["cell"] = counts["cell"].astype(str)
counts["tot_count"] = counts["w"] + counts["c"]

# Same cell threshold as the R script: the smallest library among the cells
# that passed the mosaic count QC (pass1 == 1 in info_raw).
info = pd.read_csv(snakemake.input.info_raw, sep="\t", skiprows=13)
passed = info.loc[info["pass1"] == 1, "good"]
if passed.empty:
    raise ValueError(f"no cell with pass1 == 1 in {snakemake.input.info_raw}")
min_reads = int(passed.min()) - 1

norm = library_normalization(counts, min_reads=min_reads, plot=False)
print(
    f"min_reads={min_reads}: {norm['cell'].nunique()} / {counts['cell'].nunique()} cells kept"
)

cols = ["chrom", "start", "end", "sample", "cell", "w", "c", "tot_count", "class"]
norm[cols].to_csv(snakemake.output.counts_scaled, index=False)
