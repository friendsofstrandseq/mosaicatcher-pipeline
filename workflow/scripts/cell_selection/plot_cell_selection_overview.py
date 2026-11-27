"""Multi-plate cell selection overview of a data_location.

Input:  plots/plate/cell_selection_well_table.tsv of every sample
        (plot_plate_cell_selection.py)
Output: plots/cell_selection_overview.pdf and .tsv (one row per sample)

Page 1: summary table (cells, ashleys-qc, EVAs pass/borderline/fail, union,
intersection, decision, cells kept, cells used by MosaiCatcher, median reads).
Next pages: one row of mini-plates per sample (ashleys-qc probability, EVAs
status, reads, final selection), 6 samples per page.
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import ListedColormap, LogNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, snakemake.scriptdir)
import plate_layout as pl  # noqa: E402

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

tables = [pd.read_csv(f, sep="\t") for f in snakemake.input.well_tables]
tables = [t for t in tables if len(t)]
summary = pd.DataFrame([pl.summarise(t) for t in tables])
if len(summary):
    summary = summary.sort_values("sample")
summary.to_csv(snakemake.output.tsv, sep="\t", index=False)

# Common colour scale for reads and EVAs categories across samples
all_cells = pd.concat(tables) if tables else pd.DataFrame()
if "evas_status" in all_cells.columns:
    all_cells["evas_category"] = all_cells.apply(pl.evas_category, axis=1)
    evas_colors = pl.category_colors(
        pl.sort_categories(all_cells["evas_category"].dropna().unique())
    )
else:
    evas_colors = {}
reads = (
    all_cells["mc_mapped"]
    if "mc_mapped" in all_cells.columns
    else pd.Series(dtype=float)
)
reads = reads[reads > 0]
reads_norm = LogNorm(vmin=reads.min(), vmax=reads.max()) if len(reads) else None

PANELS = ["ashleys-qc probability", "EVAs status", "reads (mapped)", "final selection"]
PER_PAGE = 6


def fmt(v, digits=0):
    return "-" if pd.isna(v) else (f"{v:,.{digits}f}" if digits else f"{int(v):,}")


def mini(ax, t, panel, n_rows, n_cols):
    placed = t.loc[t["row"].notna()]
    if placed.empty:
        return pl.not_run(
            ax, "", n_rows, n_cols, fontsize=6, ticks=False, text=pl.UNPLACED
        )
    if panel == "ashleys-qc probability":
        if "ashleys_probability" not in t.columns:
            return pl.not_run(ax, "", n_rows, n_cols, fontsize=6, ticks=False)
        cm = matplotlib.colormaps["RdYlGn"].copy()
        cm.set_bad(pl.NA_COLOR)
        ax.imshow(
            np.ma.masked_invalid(
                pl.grid_of(placed, "ashleys_probability", n_rows, n_cols)
            ),
            cmap=cm,
            vmin=0,
            vmax=1,
        )
    elif panel == "EVAs status":
        if "evas_status" not in t.columns:
            return pl.not_run(ax, "", n_rows, n_cols, fontsize=6, ticks=False)
        cats = list(evas_colors)
        idx = {c: i for i, c in enumerate(cats)}
        placed = placed.assign(_v=placed.apply(pl.evas_category, axis=1).map(idx))
        cmap = ListedColormap([evas_colors[c] for c in cats])
        cmap.set_bad(pl.NA_COLOR)
        ax.imshow(
            np.ma.masked_invalid(pl.grid_of(placed, "_v", n_rows, n_cols)),
            cmap=cmap,
            vmin=-0.5,
            vmax=len(cats) - 0.5,
        )
    elif panel == "reads (mapped)":
        if "mc_mapped" not in t.columns or reads_norm is None:
            return pl.not_run(ax, "", n_rows, n_cols, fontsize=6, ticks=False)
        g = pl.grid_of(placed, "mc_mapped", n_rows, n_cols)
        cm = matplotlib.colormaps["viridis"].copy()
        cm.set_bad(pl.NA_COLOR)
        ax.imshow(
            np.ma.masked_invalid(np.where(g > 0, g, np.nan)), cmap=cm, norm=reads_norm
        )
    else:
        cmap = ListedColormap(pl.KEEP_COLORS)
        cmap.set_bad(pl.NA_COLOR)
        ax.imshow(
            np.ma.masked_invalid(pl.grid_of(placed, "prediction", n_rows, n_cols)),
            cmap=cmap,
            vmin=0,
            vmax=1,
        )
    pl.draw_grid(ax, "", n_rows, n_cols, ticks=False)


with PdfPages(snakemake.output.pdf) as pdf:
    # Page 1: summary table
    cols = [
        ("sample", "sample", None),
        ("cells", "cells", 0),
        ("ashleys", "ashleys_kept", 0),
        ("EVAs pass", "evas_pass", 0),
        ("border", "evas_borderline", 0),
        ("fail", "evas_fail", 0),
        ("union", "union", 0),
        ("inter.", "intersection", 0),
        ("decision", "decision", None),
        ("kept", "kept", 0),
        ("used by MC", "used_by_mosaicatcher", 0),
        ("median reads", "median_reads", 0),
        ("median dup.", "median_dupl_rate", 2),
    ]
    rows = [
        [
            str(r[c]) if d is None else (fmt(r[c], d) if d else fmt(r[c]))
            for _, c, d in cols
        ]
        for _, r in summary.iterrows()
    ]
    height = max(3, 0.32 * (len(rows) + 3))
    fig, ax = plt.subplots(figsize=(16, height))
    ax.axis("off")
    if rows:
        t = ax.table(
            cellText=rows,
            colLabels=[c[0] for c in cols],
            loc="upper center",
            cellLoc="center",
        )
        t.auto_set_font_size(False)
        t.set_fontsize(7)
        t.scale(1, 1.3)
        t.auto_set_column_width(list(range(len(cols))))
    ax.set_title(
        f"Cell selection overview: {len(summary)} samples, {int(summary['cells'].sum()) if len(summary) else 0} cells",
        fontsize=11,
    )
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Mini-plates, PER_PAGE samples per page
    ordered = sorted(tables, key=lambda t: str(t["sample"].iloc[0]))
    for start in range(0, len(ordered), PER_PAGE):
        chunk = ordered[start : start + PER_PAGE]
        fig, axes = plt.subplots(
            len(chunk), len(PANELS), figsize=(14, 2.6 * len(chunk) + 1), squeeze=False
        )
        for i in range(len(chunk)):
            for j in range(len(PANELS)):
                axes[i, j].axis("off")
        for i, t in enumerate(chunk):
            t, plate, n_rows, n_cols = pl.place_wells(t)
            s = pl.summarise(t)
            for j, panel in enumerate(PANELS):
                ax = axes[i, j]
                ax.axis("on")
                mini(ax, t, panel, n_rows, n_cols)
                if i == 0:
                    ax.set_title(panel, fontsize=9)
            n_placed = int(t["row"].notna().sum())
            label = (
                f"{s['sample']}\n{s['cells']} cells ({plate}-well"
                + (
                    f", {s['cells'] - n_placed} not placed"
                    if n_placed < s["cells"]
                    else ""
                )
                + ")\n"
                f"kept {s['kept']} ({s['decision']})\n"
                f"ashleys {fmt(s['ashleys_kept'])} | EVAs {fmt(s['evas_pass'])}"
            )
            axes[i, 0].set_ylabel(
                label, fontsize=6.5, rotation=0, ha="right", va="center", labelpad=8
            )
        handles = [Patch(color=c, label=k) for k, c in evas_colors.items()]
        handles += [
            Patch(color=pl.KEEP_COLORS[1], label="kept"),
            Patch(color=pl.KEEP_COLORS[0], label="removed"),
        ]
        fig.legend(
            handles=handles,
            loc="lower center",
            ncol=min(5, len(handles)),
            fontsize=7,
            frameon=False,
        )
        fig.suptitle(
            "Cell selection per plate (ashleys-qc: red-green 0-1; reads: log scale, common to all samples)",
            fontsize=10,
        )
        fig.tight_layout(
            rect=(
                0.12,
                0.6 / (2.6 * len(chunk) + 1),
                1,
                1 - 0.4 / (2.6 * len(chunk) + 1),
            )
        )
        pdf.savefig(fig)
        plt.close(fig)

print(summary.to_string(index=False))
