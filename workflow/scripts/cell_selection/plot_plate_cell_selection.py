"""Plate report of one sample: cell selection (ashleys-qc, EVAs, decision) and
MosaiCatcher per-cell metrics.

Inputs:  cell_selection/labels_detailed.tsv (decide_cell_selection.py),
         counts/{sample}.info_raw (mosaic_count), the labels.tsv used by
         MosaiCatcher, and the SV calls if they already exist (optional).
Outputs: plots/plate/cell_selection_plate.pdf and a per-well table.

Pages: 1 = ashleys-qc probability | EVAs status | final decision;
       2 = MosaiCatcher layer: reads, duplicate rate, coverage pass (pass1),
           selection used by MosaiCatcher, EVAs noisy fraction, SV calls/cell;
       3 = ashleys-qc vs EVAs (only when both ran).
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

sample = snakemake.wildcards.sample
df = pd.read_csv(snakemake.input.detailed, sep="\t")
df["sample"] = sample
df["cell_key"] = df["cell"].map(pl.cell_key)
df = df.merge(pl.read_info(snakemake.input.info), on="cell_key", how="left")
df = df.merge(pl.read_used_labels(snakemake.input.labels), on="cell_key", how="left")
sv = pl.read_sv_counts(snakemake.params.sv_calls)
if sv is not None:
    df = df.merge(sv, on="cell_key", how="left")
    df["mc_sv_calls"] = df["mc_sv_calls"].fillna(0)
decision = df["decision"].iloc[0] if len(df) else "none"

df, plate, n_rows, n_cols = pl.place_wells(df)
valid = df["row"].notna()
placed = df.loc[valid]
has_ashleys = "ashleys_probability" in df.columns
has_evas = "evas_status" in df.columns
if has_evas:
    df["evas_category"] = df.apply(pl.evas_category, axis=1)
    placed = df.loc[valid]
small = plate == 384  # smaller text in 384-well plates
text_size = 3 if small else 4.5


def categorical(ax, column, colors, title, legend_counts=True):
    cats = list(colors)
    index = {c: i for i, c in enumerate(cats)}
    grid = np.full((n_rows, n_cols), np.nan)
    for r, c, v in zip(placed["row"], placed["col"], placed[column]):
        if v in index:
            grid[int(r), int(c)] = index[v]
    cmap = ListedColormap([colors[c] for c in cats])
    cmap.set_bad(pl.NA_COLOR)
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=-0.5, vmax=len(cats) - 0.5)
    counts = df[column].value_counts()
    ax.legend(
        handles=[
            Patch(
                color=colors[c],
                label=f"{c} ({counts.get(c, 0)})" if legend_counts else c,
            )
            for c in cats
        ],
        fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
    )
    pl.draw_grid(ax, title, n_rows, n_cols)


def continuous(
    fig, ax, column, title, cmap="viridis", log=False, fmt=None, vmin=None, vmax=None
):
    if column not in df.columns or df[column].isna().all():
        pl.not_run(ax, title + " (not available)", n_rows, n_cols)
        return
    grid = pl.grid_of(placed, column, n_rows, n_cols)
    cm = matplotlib.colormaps[cmap].copy()
    cm.set_bad(pl.NA_COLOR)
    if log:
        grid = np.where(grid > 0, grid, np.nan)
        norm = LogNorm(
            vmin=np.nanmin(grid) if np.isfinite(grid).any() else 1,
            vmax=np.nanmax(grid) if np.isfinite(grid).any() else 10,
        )
        im = ax.imshow(np.ma.masked_invalid(grid), cmap=cm, norm=norm)
    else:
        im = ax.imshow(np.ma.masked_invalid(grid), cmap=cm, vmin=vmin, vmax=vmax)
    if fmt and not small:
        for r, c, v in zip(placed["row"], placed["col"], placed[column]):
            if pd.notna(v):
                ax.text(c, r, fmt(v), ha="center", va="center", fontsize=text_size)
    pl.draw_grid(ax, title, n_rows, n_cols)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)


def kept_removed(ax, column, title, marks=None):
    grid = pl.grid_of(placed, column, n_rows, n_cols)
    cmap = ListedColormap(pl.KEEP_COLORS)
    cmap.set_bad(pl.NA_COLOR)
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=0, vmax=1)
    handles = [
        Patch(color=pl.KEEP_COLORS[1], label="kept"),
        Patch(color=pl.KEEP_COLORS[0], label="removed"),
    ]
    if marks is not None:
        for _, r in placed.loc[marks].iterrows():
            ax.text(
                r["col"],
                r["row"],
                "x",
                ha="center",
                va="center",
                fontsize=6 if small else 8,
            )
        handles.append(Patch(color="white", label="x = tools disagree"))
    ax.legend(
        handles=handles,
        fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
    )
    pl.draw_grid(ax, title, n_rows, n_cols)


with PdfPages(snakemake.output.pdf) as pdf:
    # Page 1: cell selection
    fig, axes = plt.subplots(1, 3, figsize=(20, 6 if not small else 7))
    if has_ashleys:
        n_ok = int(df["ashleys_prediction"].sum())
        continuous(
            fig,
            axes[0],
            "ashleys_probability",
            f"ashleys-qc probability ({n_ok}/{len(df)} >= threshold)",
            cmap="RdYlGn",
            vmin=0,
            vmax=1,
            fmt=lambda v: f"{v:.2f}",
        )
    else:
        pl.not_run(axes[0], "ashleys-qc probability", n_rows, n_cols)
    if has_evas:
        colors = pl.category_colors(
            pl.sort_categories(df["evas_category"].dropna().unique())
        )
        categorical(
            axes[1],
            "evas_category",
            colors,
            f"EVAs status ({int(df['evas_selected'].sum())}/{len(df)} selected)"
            + (
                "\nborderline: K keep / R review / D discard (recommendation)"
                if "borderline_recommendation" in df.columns
                else ""
            ),
        )
        if "borderline_recommendation" in df.columns:
            for _, r in placed.loc[placed["evas_status"] == "borderline"].iterrows():
                axes[1].text(
                    r["col"],
                    r["row"],
                    pl.BORDER_MARK.get(r["borderline_recommendation"], "?"),
                    ha="center",
                    va="center",
                    fontsize=5 if small else 7,
                    fontweight="bold",
                )
    else:
        pl.not_run(axes[1], "EVAs status", n_rows, n_cols)
    disagree = (
        (placed["ashleys_prediction"] != placed["evas_selected"])
        if has_ashleys and has_evas
        else None
    )
    kept_removed(
        axes[2],
        "prediction",
        f"Final selection, decision={decision} ({int(df['prediction'].sum())}/{len(df)} kept)",
        marks=disagree.index[disagree] if disagree is not None else None,
    )
    unplaced = len(df) - int(valid.sum())
    note = (
        f" | {unplaced} cells not placed: {pl.UNPLACED.replace(chr(10), ' ')}"
        if unplaced
        else ""
    )
    fig.suptitle(f"{sample} | cell selection | {plate}-well plate{note}", fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Page 2: MosaiCatcher layer
    fig, axes = plt.subplots(2, 3, figsize=(20, 10 if not small else 12))
    reads = df["mc_mapped"] if "mc_mapped" in df.columns else pd.Series(dtype=float)
    continuous(
        fig,
        axes[0, 0],
        "mc_mapped",
        f"Reads (mapped, median {reads.median():,.0f})",
        cmap="viridis",
        log=True,
    )
    continuous(
        fig,
        axes[0, 1],
        "mc_dupl_rate",
        (
            f"Duplicate rate (median {df['mc_dupl_rate'].median():.2f})"
            if "mc_dupl_rate" in df.columns
            else "Duplicate rate"
        ),
        cmap="magma_r",
        vmin=0,
        vmax=1,
        fmt=lambda v: f"{v:.2f}",
    )
    if "mc_pass1" in df.columns:
        df["mc_pass1_label"] = df["mc_pass1"].map(
            {1: "enough coverage", 0: "low coverage"}
        )
        placed = df.loc[valid]
        categorical(
            axes[0, 2],
            "mc_pass1_label",
            {"enough coverage": "#1a9850", "low coverage": "#bdbdbd"},
            "mosaicatcher count: coverage (pass1)",
        )
    else:
        pl.not_run(axes[0, 2], "coverage (pass1)", n_rows, n_cols)
    if "mc_used" in df.columns and df["mc_used"].notna().any():
        kept_removed(
            axes[1, 0],
            "mc_used",
            f"Selection used by MosaiCatcher ({int(df['mc_used'].sum())}/{len(df)})",
        )
    else:
        pl.not_run(axes[1, 0], "Selection used by MosaiCatcher", n_rows, n_cols)
    if has_evas and "evas_noisy_fraction" in df.columns:
        continuous(
            fig,
            axes[1, 1],
            "evas_noisy_fraction",
            "EVAs noisy fraction",
            cmap="cividis",
        )
    else:
        pl.not_run(axes[1, 1], "EVAs noisy fraction", n_rows, n_cols)
    if "mc_sv_calls" in df.columns:
        continuous(
            fig,
            axes[1, 2],
            "mc_sv_calls",
            "SV calls per cell (stringent)",
            cmap="Blues",
            fmt=lambda v: f"{v:.0f}",
        )
    else:
        pl.not_run(axes[1, 2], "SV calls per cell (no SV calls yet)", n_rows, n_cols)
    fig.suptitle(
        f"{sample} | MosaiCatcher per-cell metrics | {plate}-well plate", fontsize=11
    )
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Page 3: tool comparison (only when both ran)
    if has_ashleys and has_evas:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 5))
        order = list(colors)
        rng = np.random.default_rng(0)
        for i, c in enumerate(order):
            y = df.loc[df["evas_category"] == c, "ashleys_probability"]
            ax1.scatter(
                rng.normal(i, 0.08, len(y)), y, s=12, color=colors[c], alpha=0.8
            )
        ax1.axhline(
            float(snakemake.params.ashleys_threshold), color="grey", ls="--", lw=0.8
        )
        ax1.set_xticks(range(len(order)))
        ax1.set_xticklabels(order, rotation=30, ha="right", fontsize=7)
        ax1.set_ylabel("ashleys-qc probability")
        ax1.set_title("ashleys-qc probability by EVAs status", fontsize=9)
        if "mc_mapped" in df.columns:
            for c in order:
                sel = df["evas_category"] == c
                ax2.scatter(
                    df.loc[sel, "mc_mapped"].clip(lower=1),
                    df.loc[sel, "ashleys_probability"],
                    s=12,
                    color=colors[c],
                    alpha=0.8,
                    label=c,
                )
            ax2.set_xscale("log")
            ax2.set_xlabel("reads (mapped)")
            ax2.set_ylabel("ashleys-qc probability")
            ax2.legend(fontsize=6, frameon=False)
            ax2.set_title(
                "reads vs ashleys-qc probability, colour = EVAs status", fontsize=9
            )
        else:
            ax2.axis("off")
        tab = pd.crosstab(
            df["ashleys_prediction"].map({1: "ashleys keep", 0: "ashleys remove"}),
            df["evas_status"],
            margins=True,
        )
        ax3.axis("off")
        t = ax3.table(
            cellText=tab.values,
            rowLabels=tab.index,
            colLabels=tab.columns,
            loc="center",
            cellLoc="center",
        )
        t.auto_set_font_size(False)
        t.set_fontsize(8)
        t.scale(1, 1.6)
        ax3.set_title("ashleys-qc vs EVAs (cells)", fontsize=9)
        fig.suptitle(f"{sample} | ashleys-qc vs EVAs", fontsize=11)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    # Page 4: EVAs borderline cells and their recommendation, to review
    border = pl.borderline_rows(df)
    if len(border):
        header = [
            "well",
            "cell",
            "ashleys",
            "noisy fr.",
            "reads",
            "recommendation",
            "rationale",
        ]
        rows = [
            [
                r["well"],
                pl.cell_key(r["cell"]),
                f"{r['ashleys_probability']:.2f}" if has_ashleys else "-",
                f"{r['evas_noisy_fraction']:.3f}",
                f"{r['evas_tot_count']:,.0f}",
                r["borderline_recommendation"],
                str(r["borderline_rationale"]).replace(" | ", "\n", 1)[:300],
            ]
            for _, r in border.iterrows()
        ]
        pl.table_pages(
            pdf,
            plt,
            header,
            rows,
            f"{sample} | EVAs borderline cells ({len(rows)}): recommendation "
            f"(applied only with evas_borderline: recommend)",
            colors=[pl.BORDER_COLORS.get(r[5]) for r in rows],
            per_page=12,
            widths=[0.04, 0.17, 0.05, 0.06, 0.07, 0.08, 0.53],
        )

cols = ["cell", "sample", "well", "well_number", "prediction", "decision"]
cols += [
    c
    for c in df.columns
    if c.startswith(("ashleys_", "evas_", "mc_", "borderline_"))
    and c not in ("mc_pass1_label", "evas_category")
]
df[cols].sort_values("well_number").to_csv(
    snakemake.output.well_table, sep="\t", index=False
)
print(f"plate={plate}, cells={len(df)}, placed={int(valid.sum())}")
