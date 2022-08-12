import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# SETTINGS
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 2.50


# LOADING DATAFRAME
df = pd.read_csv(snakemake.input.ploidy_detailled, sep="\t")
df = df.loc[df["#chrom"] != "genome"]

# GETTING CHR LIST & ORDER CATEGORICALLY
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX"]
df["#chrom"] = pd.Categorical(df["#chrom"], categories=chroms, ordered=True)
df = df.sort_values(by=["#chrom", "start"])

# STARTING SUBPLOTS
f, ax = plt.subplots(ncols=len(df["#chrom"].unique()), figsize=(3 * len(df["#chrom"].unique()), 35))

# ITERATE OVER CHROM
for i, chrom in enumerate(df["#chrom"].unique().tolist()):
    # PLOTTING MAIN FIGURE
    ax[i].plot(df.loc[df["#chrom"] == chrom].ploidy_estimate, df.loc[df["#chrom"] == chrom].start, lw=4, color="black")

    # CUSTOMISATION
    ax[i].set_xlabel("{}".format(chrom), fontsize=30)
    ax[i].set_ylim(0, df.start.max())
    ax[i].set_xlim(0, 6)

    # ADDING VERTICAL RED DASHED LINE
    ax[i].axvline(2, ymax=df.loc[df["#chrom"] == chrom].start.max() / df.start.max(), ls="--", lw=2, color="red")

    # GRID CUSTOMISATION
    ax[i].tick_params(axis="x", which="major", labelsize=20)
    major_ticks = np.arange(0, 7, 1)
    ax[i].set_xticks(major_ticks)
    ax[i].grid(which="both")
    for axe in ["top", "bottom", "left", "right"]:
        ax[i].spines[axe].set_linewidth(2)
        ax[i].spines[axe].set_color("black")
    ax[i].grid(axis="both", which="major")
    if i == 0:
        ax[i].tick_params(axis="y", which="major", labelsize=20)
        ax[i].set_ylabel("Position (Mbp)", fontsize=30)
    else:
        ax[i].get_yaxis().set_visible(False)
f.suptitle("Sample: {}".format("RPE-BM510"), fontsize=50)
plt.tight_layout()
plt.subplots_adjust(top=0.95)
plt.savefig(snakemake.output[0], dpi=300)
