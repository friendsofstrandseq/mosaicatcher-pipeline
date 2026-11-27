"""Combine per-cell QC predictions (ashleys-qc and/or EVAs) into labels.tsv.

snakemake.input:
    info     mosaicatcher count info (counts/{sample}.info_raw): every cell
    ashleys  (optional) ashleys-qc labels after the positive/negative control
             correction: cell, prediction, probability, sample
    evas     (optional) EVAs predictions: cell, prediction (1 / 0.5 / 0), reason, ...

snakemake.output:
    labels    same format as before (cell, prediction, probability, sample); the
              rest of the pipeline keeps the cells with prediction == 1
    detailed  one row per cell: final prediction, decision mode, the columns
              of each tool (ashleys_*, evas_*) and, for EVAs borderline cells,
              borderline_recommendation (keep / review / discard) and
              borderline_rationale (borderline.py)
"""

import sys

import pandas as pd

sys.path.insert(0, snakemake.scriptdir)
import borderline  # noqa: E402

SUFFIX = ".sort.mdup.bam"

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sample = snakemake.wildcards.sample
decision = snakemake.params.decision
threshold = float(snakemake.params.ashleys_threshold)
borderline_mode = snakemake.params.evas_borderline  # exclude | recommend | include
inputs = dict(snakemake.input.items())


def with_suffix(cell):
    cell = str(cell)
    return cell if cell.endswith(SUFFIX) else cell + SUFFIX


info = pd.read_csv(inputs["info"], sep="\t", skiprows=13)
cells = sorted(info["cell"].map(with_suffix).unique())
print(f"{sample}: {len(cells)} cells in {inputs['info']}")

# Base table. With ashleys-qc, start from its labels so that labels.tsv keeps
# exactly the rows and columns produced before this step existed.
if "ashleys" in inputs:
    labels = pd.read_csv(inputs["ashleys"], sep="\t")
    labels["cell"] = labels["cell"].map(with_suffix)
    extra = sorted(set(cells) - set(labels["cell"]))
    if extra:
        print(f"{len(extra)} cells missing from the ashleys-qc labels, added as 0")
        labels = pd.concat(
            [
                labels,
                pd.DataFrame(
                    {
                        "cell": extra,
                        "prediction": 0,
                        "probability": 0.0,
                        "sample": sample,
                    }
                ),
            ],
            ignore_index=True,
        )
else:
    labels = pd.DataFrame(
        {"cell": cells, "prediction": 0, "probability": 0.0, "sample": sample}
    )

detail = labels[["cell"]].copy()
detail["sample"] = sample

if "ashleys" in inputs:
    detail["ashleys_probability"] = labels["probability"].astype(float).values
    # Same rule as the historical tune_predictions_based_on_threshold step
    detail["ashleys_prediction"] = (detail["ashleys_probability"] >= threshold).astype(
        int
    )
    print(
        f"ashleys-qc: {int(detail['ashleys_prediction'].sum())}/{len(detail)} cells "
        f"with probability >= {threshold}"
    )

if "evas" in inputs:
    ev = pd.read_csv(inputs["evas"], sep="\t")
    ev["cell"] = ev["cell"].map(with_suffix)
    cols = [
        c
        for c in [
            "prediction",
            "reason",
            "noisy_fraction",
            "mean_interpoint_distance",
            "n_extreme_wf",
            "no_of_bins_below_20_percent",
            "tot_count",
        ]
        if c in ev.columns
    ]
    ev = ev[["cell"] + cols].rename(columns={c: f"evas_{c}" for c in cols})
    detail = detail.merge(ev, on="cell", how="left")
    status = {1.0: "pass", 0.5: "borderline", 0.0: "fail"}
    detail["evas_status"] = (
        detail["evas_prediction"].astype(float).map(status).fillna("missing")
    )
    if "evas_reason" not in detail.columns:
        detail["evas_reason"] = ""
    detail["evas_reason"] = detail["evas_reason"].fillna("")
    detail.loc[detail["evas_status"] == "missing", "evas_reason"] = "not in EVAs output"
    mc = info.assign(cell=info["cell"].map(with_suffix))
    rec, why = borderline.recommend(
        detail,
        mc,
        ashleys_threshold=threshold,
        extra_args=snakemake.params.evas_extra_args,
    )
    detail["borderline_recommendation"] = rec.values
    detail["borderline_rationale"] = why.values
    selected = detail["evas_status"] == "pass"
    is_border = detail["evas_status"] == "borderline"
    if borderline_mode == "include":
        selected |= is_border
    elif borderline_mode == "recommend":
        selected |= is_border & (detail["borderline_recommendation"] == "keep")
    detail["evas_selected"] = selected.astype(int)
    detail["evas_borderline_mode"] = borderline_mode
    print(
        "EVAs: "
        + ", ".join(f"{k}={v}" for k, v in detail["evas_status"].value_counts().items())
        + f"; borderline recommendations: "
        + ", ".join(
            f"{k}={v}"
            for k, v in detail.loc[is_border, "borderline_recommendation"]
            .value_counts()
            .items()
        )
        + f"; evas_borderline={borderline_mode}"
    )

if decision == "ashleys":
    keep = detail["ashleys_prediction"]
elif decision == "evas":
    keep = detail["evas_selected"]
elif decision == "union":
    keep = detail["ashleys_prediction"] | detail["evas_selected"]
elif decision == "intersection":
    keep = detail["ashleys_prediction"] & detail["evas_selected"]
elif decision == "none":
    # Predictions are computed and reported, every cell is kept
    keep = pd.Series(1, index=detail.index)
else:
    raise ValueError(f"Unknown cell_selection decision: {decision}")

detail["prediction"] = keep.astype(int).values
detail["decision"] = decision
labels["prediction"] = detail["prediction"].values
if "ashleys" not in inputs:
    # No ashleys-qc probability: report the EVAs score (1 / 0.5 / 0), or 1.0
    labels["probability"] = (
        detail["evas_prediction"].astype(float).fillna(0.0).values
        if "evas" in inputs
        else 1.0
    )

print(
    f"decision={decision}: {int(labels['prediction'].sum())}/{len(labels)} cells kept"
)

labels.sort_values(by="cell").to_csv(snakemake.output.labels, sep="\t", index=False)
first = ["cell", "sample", "prediction", "decision"]
detail[first + [c for c in detail.columns if c not in first]].sort_values(
    by="cell"
).to_csv(snakemake.output.detailed, sep="\t", index=False)
