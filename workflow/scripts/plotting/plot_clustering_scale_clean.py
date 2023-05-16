import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

log = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log

# Categorical mapping
# d = {
#     "none": 0,
#     "del_h1": 1,
#     "del_h2": 2,
#     "del_hom": 3,
#     "dup_h1": 4,
#     "dup_h2": 5,
#     "dup_hom": 6,
#     "inv_h1": 7,
#     "inv_h2": 8,
#     "inv_hom": 9,
#     "idup_h1": 10,
#     "idup_h2": 11,
#     "complex": 12,
# }

d = {
    "none": 1,
    "del_h1": 2,
    "del_h2": 3,
    "del_hom": 4,
    "dup_h1": 5,
    "dup_h2": 6,
    "dup_hom": 7,
    "inv_h1": 8,
    "inv_h2": 9,
    "inv_hom": 10,
    "idup_h1": 11,
    "idup_h2": 12,
    "complex": 13,
}

# Colors
colors = {
    "none": "#F8F8F8",
    "del_h1": "#77AADD",
    "del_h2": "#4477AA",
    "del_hom": "#114477",
    "dup_h1": "#CC99BB",
    "dup_h2": "#AA4488",
    "dup_hom": "#771155",
    "inv_h1": "#DDDD77",
    "inv_h2": "#AAAA44",
    "inv_hom": "#777711",
    "idup_h1": "#DDAA77",
    "idup_h2": "#AA7744",
    "complex": "#774411",
}

# Read SV file
df = pd.read_csv(snakemake.input.sv_calls, sep="\t")
# df = pd.read_csv("../stringent_filterTRUE.tsv", sep="\t")
df["ID"] = df["chrom"] + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)

# Read 200kb bins file
binbed = pd.read_csv(
    # "../bin_200kb_all.bed",
    snakemake.input.binbed,
    sep="\t",
    names=["chrom", "start", "end", "bin_id"],
)
binbed["ID"] = binbed["chrom"] + "_" + binbed["start"].astype(str) + "_" + binbed["end"].astype(str)

# Turn chrom into categorical
binbed["chrom"] = pd.Categorical(
    binbed["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)

# Sort & filter out chrY #TMP / can be changed
binbed = binbed.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)
# binbed = binbed.loc[~binbed["chrom"].isin(["chrY"])]

# Instanciate final list
l = list()


def process_row(r):
    """Get all bins from binbed that overlap SV call

    Args:
        r (pandas row)
    """
    tmp_r = binbed.loc[(binbed["chrom"] == r["chrom"]) & (binbed["start"] >= r["start"]) & (binbed["end"] <= r["end"])]
    tmp_r["cell"] = r["cell"]
    tmp_r["sv_call_name"] = r["sv_call_name"]
    tmp_r["af"] = r["af"]
    tmp_r["llr_to_ref"] = r["llr_to_ref"]
    # Append result to list
    l.append(tmp_r)


# Apply & loop on each temporary cell dataframe created
def process_sv(tmp_df):
    tmp_df.apply(lambda r: process_row(r), axis=1)


# Create a nested pandas apply
df.groupby("cell").apply(lambda r: process_sv(r))

# Concat results
processed_df = pd.concat(l)
processed_df["ID"] = processed_df["chrom"].astype(str) + "_" + processed_df["start"].astype(str) + "_" + processed_df["end"].astype(str)

# Extract only empty bins (outer join) from binbed
binbed_not_used = binbed.loc[~binbed["ID"].isin(processed_df.ID.unique().tolist())]
# Concat with previously created dataframe
concat_df = pd.concat([processed_df, binbed_not_used])

# Replace llr inf values by max values
concat_df.loc[concat_df["llr_to_ref"] == np.inf, "llr_to_ref"] = concat_df.loc[concat_df["llr_to_ref"] != np.inf]["llr_to_ref"].max()

# Pivot into matrix
pivot_concat_df = concat_df.pivot(index="ID", values="llr_to_ref", columns="cell")

# Create chrom, start, end columns in a tmp df
tmp = pivot_concat_df.reset_index().ID.str.split("_", expand=True)
tmp.columns = ["chrom", "start", "end"]
tmp["start"] = tmp["start"].astype(int)
tmp["end"] = tmp["end"].astype(int)

# Concat dfs, remove first column, sort, index
pivot_concat_df = (
    pd.concat([pivot_concat_df.reset_index(), tmp], axis=1)
    .drop(pivot_concat_df.columns[0], axis=1)
    .sort_values(by=["chrom", "start", "end"])
    .reset_index(drop=True)
)


# Read clustering index file produced from previous clustering using R ComplexHeatmap
# clustering_index_df = pd.read_csv("test.tsv", sep="\t")
clustering_index_df = pd.read_csv(snakemake.input.cluster_order_df, sep="\t")


## LLR

# Pivot df subset specific to llr
pivot_concat_df = concat_df.pivot(index="ID", values="llr_to_ref", columns="cell")
tmp = pivot_concat_df.reset_index().ID.str.split("_", expand=True)
tmp.columns = ["chrom", "start", "end"]
tmp["start"] = tmp["start"].astype(int)
tmp["end"] = tmp["end"].astype(int)
pivot_concat_df = pd.concat([pivot_concat_df.reset_index(), tmp], axis=1).drop(pivot_concat_df.columns[0], axis=1)

pivot_concat_df["chrom"] = pd.Categorical(
    pivot_concat_df["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)
pivot_concat_df = pivot_concat_df.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)

chroms = ["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"]
# chroms = chroms[:2]
# chroms = ["chr10", "chr13", "chr22"]

# Extract widths using binbed max values to specify subplots widths scaled according chrom sizes
widths = binbed.loc[binbed["chrom"].isin(chroms)].groupby("chrom")["end"].max().dropna().tolist()


# pdf = PdfPages("multipage_pdf2.pdf")
pdf = PdfPages(snakemake.output.pdf)

# Create subplots
f, axs = plt.subplots(ncols=len(chroms), figsize=(40, 20), gridspec_kw={"width_ratios": widths})

print("LLR plot")
# Iterate over chroms
for j, (chrom, ax) in enumerate(zip(chroms, axs)):
    print(chrom)
    cbar = False

    # If not chr1 = remove y axis
    if j != 0:
        ax.get_yaxis().set_visible(False)
        ax.yaxis.set_ticks_position("none")

    # If last chrom, enable cbar plot
    if j == len(chroms) - 1:
        cbar = True

    # Subset chrom data, set_index, transpose & replace NaN by 0
    data_heatmap = (
        pivot_concat_df.loc[pivot_concat_df["chrom"] == chrom].drop(["chrom", "start", "end"], axis=1).set_index("ID").T.fillna(0)
    )

    # Reorder rows based on clustering index
    data_heatmap = data_heatmap.loc[clustering_index_df.cell.values.tolist()]

    # Plot
    sns.heatmap(
        data=data_heatmap,
        ax=ax,
        vmin=0,
        vmax=concat_df.llr_to_ref.max(),
        cmap="Reds",
        cbar=cbar,
    )
    ax.xaxis.set_ticks_position("none")
    ax.set_xlabel("{}".format(chrom), fontsize=12, rotation=90)
    ax.set_xticklabels([])

plt.suptitle(
    f"Chromosome size scaled LLR heatmap (Sample : {snakemake.wildcards.sample}, Methods used: {snakemake.wildcards.method}, Filter used: {snakemake.wildcards.filter})",
    x=0.4,
    y=1.02,
    fontsize=18,
)

pdf.savefig(f)
plt.close()

## CATEGORICAL

# Map values to categorical names
concat_df["sv_call_name_map"] = concat_df["sv_call_name"].map(d)

# Pivot df subset specific to sv_call_name
pivot_concat_df = concat_df.pivot(index="ID", values="sv_call_name_map", columns="cell")
tmp = pivot_concat_df.reset_index().ID.str.split("_", expand=True)
tmp.columns = ["chrom", "start", "end"]
tmp["start"] = tmp["start"].astype(int)
tmp["end"] = tmp["end"].astype(int)
pivot_concat_df = pd.concat([pivot_concat_df.reset_index(), tmp], axis=1).drop(pivot_concat_df.columns[0], axis=1)

pivot_concat_df["chrom"] = pd.Categorical(
    pivot_concat_df["chrom"],
    categories=["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"],
    ordered=True,
)
pivot_concat_df = pivot_concat_df.sort_values(by=["chrom", "start", "end"]).reset_index(drop=True)


chroms = ["chr{}".format(e) for e in range(1, 23)] + ["chrX", "chrY"]
# chroms = ["chr10", "chr13"]
# chroms = chroms[:2]

widths = binbed.loc[binbed["chrom"].isin(chroms)].groupby("chrom")["end"].max().dropna().tolist()

f, axs = plt.subplots(ncols=len(chroms), figsize=(30, 15), dpi=50, gridspec_kw={"width_ratios": widths})

print("Categorical plot")
for j, (chrom, ax) in enumerate(zip(chroms, axs)):
    print(chrom)
    data_heatmap = (
        pivot_concat_df.loc[pivot_concat_df["chrom"] == chrom].drop(["chrom", "start", "end"], axis=1).set_index("ID").T.fillna(0)
    )
    data_heatmap = data_heatmap.loc[clustering_index_df.cell.values.tolist()]
    sns.heatmap(data=data_heatmap, ax=ax, vmin=0, cbar=False, cmap=list(colors.values()))
    ax.xaxis.set_ticks_position("none")

    if j != 0:
        ax.get_yaxis().set_visible(False)
        ax.yaxis.set_ticks_position("none")

    ax.set_xlabel("{}".format(chrom), fontsize=12, rotation=90)
    ax.set_xticklabels([])

custom_lines = [Line2D([0], [0], color=v, lw=12) for j, (k, v) in enumerate(colors.items())]

axs[-1].legend(
    custom_lines,
    list(colors.keys()),
    bbox_to_anchor=(1 + 0.15 * len(chroms), 0.65),
    fontsize=16,
)
plt.tight_layout(rect=[0, 0, 0.95, 1])

plt.suptitle(
    f"Chromosome size scaled categorical heatmap (Sample : {snakemake.wildcards.sample}, Methods used: {snakemake.wildcards.method}, Filter used: {snakemake.wildcards.filter})",
    x=0.4,
    y=1.02,
    fontsize=18,
)

pdf.savefig(f)
plt.close()

pdf.close()
