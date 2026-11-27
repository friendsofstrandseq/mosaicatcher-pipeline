"""GC bias correction with strandtools (replaces GC_correction.R).

LOWESS fit of log(count / median) against bin GC content on a GC-stratified
subsample, then division of each bin by the fitted trend
(strandtools.normalization.gc_correction).
"""

import sys
import types

from tqdm import std as tqdm_std

# strandtools imports tqdm.notebook, which needs ipywidgets: use the console bar.
sys.modules["tqdm.notebook"] = types.SimpleNamespace(tqdm=tqdm_std.tqdm)

import matplotlib

matplotlib.use("Agg")

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from strandtools.normalization import gc_correction

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

counts = pd.read_csv(snakemake.input.counts_scaled)
counts["cell"] = counts["cell"].astype(str)

# The pipeline GC matrices name the column "%GC"; strandtools expects "GC%".
gc = pd.read_csv(snakemake.params.gc_matrix, sep="\t")
gc_cols = [c for c in gc.columns if "GC" in c]
if not gc_cols:
    raise ValueError(f"no GC column in {snakemake.params.gc_matrix}")
gc = gc.rename(columns={gc_cols[0]: "GC%"})

# The per-GC-decile subsample is random: seed it so reruns give the same output.
np.random.seed(snakemake.params.seed)

with tempfile.TemporaryDirectory() as tmp:
    gc_path = os.path.join(tmp, "GC_matrix.tsv")
    gc.to_csv(gc_path, sep="\t", index=False)
    norm = gc_correction(
        counts,
        min_count=snakemake.params.gc_min_reads,
        n_subsample=snakemake.params.gc_n_subsample,
        plot=True,
        GC_path=gc_path,
    )

plt.gcf().savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
plt.close("all")

norm.to_csv(snakemake.output.counts_scaled_gc, index=False)
print(f"GC correction done: {len(norm)} bins, {norm['cell'].nunique()} cells")
