"""Shared helpers for the cell selection plate plots (per sample and overview).

Wells are numbered from the trailing digits of the cell name, column-major
(A1, B1, ..., H1, A2, ...), as in the ashleys-qc plate plot (plot_plate.R):
2 digits on 96-well plates (...PE20301 -> 01), 3 digits on 384-well plates.
"""

import os
import re

import numpy as np
import pandas as pd

SUFFIX = r"\.sort\.mdup\.bam$"

# EVAs status colours. Failures get one colour per reason keyword so that
# spatial patterns of a given failure mode stand out on the plate.
STATUS_COLORS = {"pass": "#2c7fb8", "borderline": "#fdae61", "missing": "#d9d9d9"}
REASON_COLORS = [  # first matching keyword wins
    ("low confidence", "#f4a582"),
    ("strand", "#8c510a"),
    ("reads", "#525252"),
    ("under", "#7b3294"),
    ("over", "#d7191c"),
]
FALLBACK = ["#e7298a", "#66a61e", "#a6761d", "#666666", "#1b9e77"]
KEEP_COLORS = ["#d7301f", "#1a9850"]  # removed, kept
NA_COLOR = "#eeeeee"


def cell_key(cell):
    """Cell name without the .sort.mdup.bam suffix (info_raw naming)."""
    return re.sub(SUFFIX, "", str(cell))


def plate_size(n_cells):
    return 384 if n_cells > 96 else 96


def place_wells(df, plate=None):
    """Add well_number, row, col, well columns; returns (df, plate, n_rows, n_cols)."""
    df = df.copy()
    plate = plate or plate_size(len(df))
    n_digits = 3 if plate == 384 else 2
    n_rows, n_cols = (16, 24) if plate == 384 else (8, 12)

    def number(cell):
        m = re.search(r"(\d{%d})$" % n_digits, cell_key(cell))
        return int(m.group(1)) if m else np.nan

    df["well_number"] = df["cell"].map(number)
    valid = df["well_number"].between(1, plate)
    df["row"] = np.where(valid, (df["well_number"] - 1) % n_rows, np.nan)
    df["col"] = np.where(valid, (df["well_number"] - 1) // n_rows, np.nan)
    df["well"] = [
        f"{chr(65 + int(r))}{int(c) + 1}" if pd.notna(r) else "Unknown"
        for r, c in zip(df["row"], df["col"])
    ]
    return df, plate, n_rows, n_cols


def evas_category(row):
    status = row.get("evas_status")
    if status != "fail":
        return status
    reason = str(row.get("evas_reason", "") or "")
    return f"fail: {reason}" if reason and reason != "nan" else "fail"


def category_colors(categories):
    fb = iter(FALLBACK * 10)
    out = {}
    for cat in categories:
        if cat in STATUS_COLORS:
            out[cat] = STATUS_COLORS[cat]
            continue
        low = str(cat).lower()
        out[cat] = next((col for key, col in REASON_COLORS if key in low), None)
        if out[cat] is None:
            out[cat] = next(fb)
    return out


def sort_categories(cats):
    return sorted(cats, key=lambda c: (c != "pass", c != "borderline", str(c)))


def read_info(path):
    """mosaicatcher count summary (info_raw): one row per cell."""
    info = pd.read_csv(path, sep="\t", comment="#")
    info["cell_key"] = info["cell"].map(cell_key)
    keep = ["cell_key", "mapped", "suppl", "dupl", "mapq", "read2", "good", "pass1"]
    info = info[[c for c in keep if c in info.columns]].copy()
    info["dupl_rate"] = info["dupl"] / info["mapped"].where(info["mapped"] > 0)
    return info.add_prefix("mc_").rename(columns={"mc_cell_key": "cell_key"})


def read_used_labels(path):
    """Selection actually used by MosaiCatcher (labels.tsv given to the SV steps)."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["cell_key", "mc_used"])
    lab = pd.read_csv(path, sep="\t")
    if "cell" not in lab.columns or "prediction" not in lab.columns:
        return pd.DataFrame(columns=["cell_key", "mc_used"])
    return pd.DataFrame(
        {"cell_key": lab["cell"].map(cell_key), "mc_used": lab["prediction"]}
    )


def read_sv_counts(path):
    """Optional: number of SV calls per cell (mosaiclassifier), if the file exists."""
    if not path or not os.path.exists(path):
        return None
    sv = pd.read_csv(path, sep="\t", usecols=lambda c: c == "cell")
    if "cell" not in sv.columns:
        return None
    counts = sv["cell"].map(cell_key).value_counts()
    return pd.DataFrame({"cell_key": counts.index, "mc_sv_calls": counts.values})


def grid_of(placed, column, n_rows, n_cols):
    grid = np.full((n_rows, n_cols), np.nan)
    for r, c, v in zip(placed["row"], placed["col"], placed[column]):
        if pd.notna(r) and pd.notna(v):
            grid[int(r), int(c)] = v
    return grid


def draw_grid(ax, title, n_rows, n_cols, fontsize=9, ticks=True):
    if ticks:
        step = 2 if n_cols > 12 else 1
        ax.set_xticks(range(0, n_cols, step))
        ax.set_xticklabels([str(i + 1) for i in range(0, n_cols, step)], fontsize=5)
        ax.set_yticks(range(0, n_rows, step))
        ax.set_yticklabels([chr(65 + i) for i in range(0, n_rows, step)], fontsize=5)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8 if n_cols <= 12 else 0.3)
    ax.tick_params(which="minor", length=0)
    ax.set_title(title, fontsize=fontsize)


def not_run(ax, title, n_rows, n_cols, fontsize=9, ticks=True, text="not run"):
    ax.imshow(np.full((n_rows, n_cols), np.nan), cmap="Greys", vmin=0, vmax=1)
    draw_grid(ax, title, n_rows, n_cols, fontsize=fontsize, ticks=ticks)
    ax.text(
        (n_cols - 1) / 2,
        (n_rows - 1) / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize + 2,
        color="grey",
    )


UNPLACED = "no well number\nin the cell names"


BORDER_MARK = {"keep": "K", "review": "R", "discard": "D"}
BORDER_COLORS = {"keep": "#1a9850", "review": "#fdae61", "discard": "#d7301f"}


def borderline_rows(df):
    """Borderline cells with their recommendation, for the report tables."""
    if "borderline_recommendation" not in df.columns:
        return df.iloc[0:0]
    b = df.loc[df["evas_status"] == "borderline"].copy()
    order = {"keep": 0, "review": 1, "discard": 2}
    return b.sort_values(
        by=["borderline_recommendation", "cell"],
        key=lambda c: c.map(order) if c.name == "borderline_recommendation" else c,
    )


def table_pages(pdf, plt, header, rows, title, colors=None, per_page=28, widths=None):
    """Draw a long table over as many pages as needed (rows = lists of str)."""
    for start in range(0, max(len(rows), 1), per_page):
        chunk = rows[start : start + per_page]
        fig, ax = plt.subplots(figsize=(16, 0.32 * (len(chunk) + 3) + 0.6))
        ax.axis("off")
        if chunk:
            t = ax.table(
                cellText=chunk,
                colLabels=header,
                loc="upper center",
                cellLoc="left",
                colWidths=widths,
            )
            t.auto_set_font_size(False)
            t.set_fontsize(6.5)
            t.scale(1, 1.25)
            if colors:
                for i, col in enumerate(colors[start : start + per_page]):
                    if col:
                        t[(i + 1, header.index("recommendation"))].set_facecolor(col)
        page = f" ({start // per_page + 1})" if len(rows) > per_page else ""
        ax.set_title(title + page, fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def summarise(df):
    """One summary row for a sample from its well table."""
    has_ashleys = "ashleys_prediction" in df.columns
    has_evas = "evas_status" in df.columns
    a = df["ashleys_prediction"] == 1 if has_ashleys else None
    e = df["evas_selected"] == 1 if has_evas else None
    row = {
        "sample": df["sample"].iloc[0] if "sample" in df.columns else "",
        "cells": len(df),
        "decision": df["decision"].iloc[0] if len(df) else "",
        "kept": int(df["prediction"].sum()),
        "ashleys_kept": int(a.sum()) if has_ashleys else np.nan,
        "evas_pass": int((df["evas_status"] == "pass").sum()) if has_evas else np.nan,
        "evas_borderline": (
            int((df["evas_status"] == "borderline").sum()) if has_evas else np.nan
        ),
        "evas_fail": int((df["evas_status"] == "fail").sum()) if has_evas else np.nan,
        **{
            f"borderline_{k}": (
                int((df["borderline_recommendation"] == k).sum())
                if "borderline_recommendation" in df.columns
                else np.nan
            )
            for k in ("keep", "review", "discard")
        },
        "union": int((a | e).sum()) if has_ashleys and has_evas else np.nan,
        "intersection": int((a & e).sum()) if has_ashleys and has_evas else np.nan,
        "used_by_mosaicatcher": (
            int(df["mc_used"].sum())
            if "mc_used" in df.columns and df["mc_used"].notna().any()
            else np.nan
        ),
        "median_reads": (
            float(df["mc_mapped"].median()) if "mc_mapped" in df.columns else np.nan
        ),
        "median_good_reads": (
            float(df["mc_good"].median()) if "mc_good" in df.columns else np.nan
        ),
        "median_dupl_rate": (
            float(df["mc_dupl_rate"].median())
            if "mc_dupl_rate" in df.columns
            else np.nan
        ),
    }
    return row
