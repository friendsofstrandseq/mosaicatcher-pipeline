"""Plate view of the cell selection: ashleys-qc probability, EVAs status and the
final decision, side by side, plus a page comparing the two tools.

Input:  cell_selection/labels_detailed.tsv (decide_cell_selection.py)
Output: plots/plate/cell_selection_plate.pdf and a per-well table
"""

import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sample = snakemake.wildcards.sample
df = pd.read_csv(snakemake.input.detailed, sep="\t")
decision = df["decision"].iloc[0] if len(df) else "none"

# Well number = trailing digits of the cell name, numbered column-major
# (A1, B1, ..., H1, A2, ...), as in the ashleys-qc plate plot (plot_plate.R).
# 96-well plates use 2 digits (...PE20301 -> 01), 384-well plates 3 digits.
plate = 384 if len(df) > 96 else 96
n_digits = 3 if plate == 384 else 2


def well_number(cell):
    name = re.sub(r"\.sort\.mdup\.bam$", "", str(cell))
    m = re.search(r"(\d{%d})$" % n_digits, name)
    return int(m.group(1)) if m else np.nan


df["well_number"] = df["cell"].map(well_number)
n_rows, n_cols = (16, 24) if plate == 384 else (8, 12)
valid = df["well_number"].between(1, plate)
df["row"] = np.nan
df["col"] = np.nan
df.loc[valid, "row"] = (df.loc[valid, "well_number"] - 1) % n_rows
df.loc[valid, "col"] = (df.loc[valid, "well_number"] - 1) // n_rows
df["well"] = [
    f"{chr(65 + int(r))}{int(c) + 1}" if pd.notna(r) else "Unknown"
    for r, c in zip(df["row"], df["col"])
]
placed = df.loc[valid].copy()

has_ashleys = "ashleys_probability" in df.columns
has_evas = "evas_status" in df.columns

# EVAs status colours. Failures get one colour per reason keyword so that
# spatial patterns of a given failure mode stand out on the plate.
STATUS_COLORS = {
    "pass": "#2c7fb8",
    "borderline": "#fdae61",
    "missing": "#d9d9d9",
}
# First matching keyword wins
REASON_COLORS = [
    ("strand", "#8c510a"),
    ("reads", "#525252"),
    ("under", "#7b3294"),
    ("over", "#d7191c"),
]
FALLBACK = ["#e7298a", "#66a61e", "#a6761d", "#666666", "#1b9e77"]


def evas_category(row):
    status = row["evas_status"]
    if status != "fail":
        return status
    reason = str(row.get("evas_reason", "") or "fail")
    return f"fail: {reason}" if reason else "fail"


def category_color(cat, fallback_iter):
    if cat in STATUS_COLORS:
        return STATUS_COLORS[cat]
    low = cat.lower()
    for key, col in REASON_COLORS:
        if key in low:
            return col
    return next(fallback_iter)


def empty_grid():
    return np.full((n_rows, n_cols), np.nan)


def draw_grid(ax, title):
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([str(i + 1) for i in range(n_cols)], fontsize=6)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([chr(65 + i) for i in range(n_rows)], fontsize=6)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)
    ax.set_title(title, fontsize=9)


def not_run(ax, title):
    ax.imshow(empty_grid(), cmap="Greys", vmin=0, vmax=1)
    draw_grid(ax, title)
    ax.text(
        (n_cols - 1) / 2,
        (n_rows - 1) / 2,
        "not run",
        ha="center",
        va="center",
        fontsize=12,
        color="grey",
    )


with PdfPages(snakemake.output.pdf) as pdf:
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.2))

    # Panel 1: ashleys probability
    ax = axes[0]
    if has_ashleys:
        grid = empty_grid()
        for _, r in placed.iterrows():
            grid[int(r["row"]), int(r["col"])] = r["ashleys_probability"]
        cmap = matplotlib.colormaps["RdYlGn"].copy()
        cmap.set_bad("#eeeeee")
        im = ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=0, vmax=1)
        for _, r in placed.iterrows():
            ax.text(
                r["col"],
                r["row"],
                f"{r['ashleys_probability']:.2f}",
                ha="center",
                va="center",
                fontsize=4.5,
            )
        n_ok = int(df["ashleys_prediction"].sum())
        draw_grid(ax, f"ashleys-qc probability ({n_ok}/{len(df)} >= threshold)")
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    else:
        not_run(ax, "ashleys-qc probability")

    # Panel 2: EVAs status / reason
    ax = axes[1]
    if has_evas:
        df["evas_category"] = df.apply(evas_category, axis=1)
        placed["evas_category"] = df.loc[valid, "evas_category"]
        cats = sorted(
            df["evas_category"].unique(),
            key=lambda c: (c != "pass", c != "borderline", c),
        )
        fb = iter(FALLBACK * 5)
        colors = [category_color(c, fb) for c in cats]
        index = {c: i for i, c in enumerate(cats)}
        grid = empty_grid()
        for _, r in placed.iterrows():
            grid[int(r["row"]), int(r["col"])] = index[r["evas_category"]]
        cmap = ListedColormap(colors)
        cmap.set_bad("#eeeeee")
        ax.imshow(
            np.ma.masked_invalid(grid), cmap=cmap, vmin=-0.5, vmax=len(cats) - 0.5
        )
        counts = df["evas_category"].value_counts()
        ax.legend(
            handles=[
                Patch(color=col, label=f"{c} ({counts.get(c, 0)})")
                for c, col in zip(cats, colors)
            ],
            fontsize=6,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            frameon=False,
        )
        draw_grid(
            ax, f"EVAs status ({int(df['evas_selected'].sum())}/{len(df)} selected)"
        )
    else:
        not_run(ax, "EVAs status")

    # Panel 3: final decision; x marks wells where the two tools disagree
    ax = axes[2]
    grid = empty_grid()
    for _, r in placed.iterrows():
        grid[int(r["row"]), int(r["col"])] = r["prediction"]
    cmap = ListedColormap(["#d7301f", "#1a9850"])
    cmap.set_bad("#eeeeee")
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=0, vmax=1)
    if has_ashleys and has_evas:
        disagree = placed["ashleys_prediction"] != placed["evas_selected"]
        for _, r in placed.loc[disagree].iterrows():
            ax.text(r["col"], r["row"], "x", ha="center", va="center", fontsize=8)
    n_keep = int(df["prediction"].sum())
    draw_grid(ax, f"Final selection, decision={decision} ({n_keep}/{len(df)} kept)")
    handles = [
        Patch(color="#1a9850", label="kept"),
        Patch(color="#d7301f", label="removed"),
    ]
    if has_ashleys and has_evas:
        handles.append(Patch(color="white", label="x = tools disagree"))
    ax.legend(
        handles=handles,
        fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
    )

    fig.suptitle(f"{sample} | cell selection | {plate}-well plate", fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    # Page 2: tool comparison (only when both ran)
    if has_ashleys and has_evas:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        order = list(cats)
        rng = np.random.default_rng(0)
        for i, c in enumerate(order):
            y = df.loc[df["evas_category"] == c, "ashleys_probability"]
            x = rng.normal(i, 0.08, len(y))
            ax1.scatter(x, y, s=12, color=colors[i], alpha=0.8)
        ax1.axhline(
            float(snakemake.params.ashleys_threshold), color="grey", ls="--", lw=0.8
        )
        ax1.set_xticks(range(len(order)))
        ax1.set_xticklabels(order, rotation=30, ha="right", fontsize=7)
        ax1.set_ylabel("ashleys-qc probability")
        ax1.set_title("ashleys-qc probability by EVAs status", fontsize=9)

        tab = pd.crosstab(
            df["ashleys_prediction"].map({1: "ashleys keep", 0: "ashleys remove"}),
            df["evas_status"],
            margins=True,
        )
        ax2.axis("off")
        t = ax2.table(
            cellText=tab.values,
            rowLabels=tab.index,
            colLabels=tab.columns,
            loc="center",
            cellLoc="center",
        )
        t.auto_set_font_size(False)
        t.set_fontsize(8)
        t.scale(1, 1.6)
        ax2.set_title("ashleys-qc vs EVAs (cells)", fontsize=9)
        fig.suptitle(f"{sample} | ashleys-qc vs EVAs", fontsize=11)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

cols = ["cell", "well", "well_number", "prediction", "decision"] + [
    c for c in df.columns if c.startswith(("ashleys_", "evas_"))
]
df[cols].sort_values("well_number").to_csv(
    snakemake.output.well_table, sep="\t", index=False
)
print(f"plate={plate}, cells={len(df)}, placed={int(valid.sum())}")
