"""Variance stabilisation with strandtools (replaces variance_stabilizing_transformation.R).

Anscombe transform for negative binomial counts. The dispersion phi is the
value on a grid that minimises the variation of the per-bin residual standard
deviation across cells (strandtools.normalization.dispersion_optimization).
"""

import sys
import types

from tqdm import std as tqdm_std

# strandtools imports tqdm.notebook, which needs ipywidgets: use the console bar.
sys.modules["tqdm.notebook"] = types.SimpleNamespace(tqdm=tqdm_std.tqdm)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from strandtools.normalization import dispersion_optimization, vst

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

counts = pd.read_csv(snakemake.input.counts_scaled_gc)
counts["cell"] = counts["cell"].astype(str)
counts["tot_count"] = counts["w"] + counts["c"]

phi = dispersion_optimization(
    counts.copy(),
    phi_min=snakemake.params.phi_min,
    phi_max=snakemake.params.phi_max,
    phi_steps=snakemake.params.phi_steps,
)
print(f"Estimated dispersion phi: {phi}")

norm = vst(counts, phi=phi, plot=True)

plt.gcf().savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
plt.close("all")

norm.to_csv(snakemake.output.counts_scaled_gc_vst, index=False)
