"""Recommendation for EVAs borderline cells: keep, review or discard.

Why a cell is borderline in EVAs (prediction rules of EVAs 0.2.1.dev, Nov 2026): it has
>= 100,000 reads, is not an over-BrdU outlier, and its noisy fraction is in
[nf_borderline_lower, nf_borderline_upper] (default [0.05, 0.08]). Below that
range a cell passes if it has >= no_strandedness_min_extreme_wf (default 100)
extreme W/C bins, above it fails ("under BrdU incorporation"). The
strandedness check is not applied to borderline cells.

The recommendation is conservative (including a non-viable cell is worse than
losing one): "keep" only when every signal agrees, "discard" when one signal is
clearly bad, "review" otherwise. Signals:

  noisy fraction    position in the borderline range: keep needs the lower
                    third (closest to pass)
  strandedness      n_extreme_wf >= the EVAs minimum: below it, the cell would
                    fail "no strandedness" even with a low noisy fraction
                    -> discard
  depth             EVAs reads (tot_count) >= 2 x 100,000
  interpoint dist.  mean_interpoint_distance inside the Tukey fences of the
                    sample's EVAs pass cells (needs >= 5 pass cells)
  ashleys-qc        probability >= max(0.8, threshold) for keep; below the
                    threshold -> discard; not run -> no keep (on 10 samples,
                    2 of the 19 cells that all the other signals would keep
                    were rejected by ashleys-qc)
  MosaiCatcher      pass1 == 1 (else discard); duplicate rate <= median of the
                    sample's EVAs pass cells + 0.10; info missing -> no keep
"""

import re

import numpy as np
import pandas as pd

DEFAULTS = {
    "nf_borderline_lower": 0.05,
    "nf_borderline_upper": 0.08,
    "no_strandedness_min_extreme_wf": 100,
    "min_reads": 100_000,
}


def evas_thresholds(extra_args=""):
    """EVAs thresholds, with the values overridden in evas_extra_args."""
    th = dict(DEFAULTS)
    for key in (
        "nf_borderline_lower",
        "nf_borderline_upper",
        "no_strandedness_min_extreme_wf",
    ):
        m = re.search(
            r"--%s[ =]([0-9.eE+-]+)" % key.replace("_", "-"), extra_args or ""
        )
        if m:
            th[key] = float(m.group(1))
    return th


def _fences(values):
    v = pd.to_numeric(values, errors="coerce").dropna()
    if len(v) < 5:
        return None
    q1, q3 = v.quantile(0.25), v.quantile(0.75)
    iqr = q3 - q1
    return q1 - 1.5 * iqr, q3 + 1.5 * iqr


def recommend(detail, info=None, ashleys_threshold=0.5, extra_args=""):
    """Return (recommendation, rationale) Series aligned on detail.index.

    detail: labels_detailed-like table (evas_* and optional ashleys_* columns)
    info:   optional MosaiCatcher info_raw table with columns cell (same naming
            as detail), mapped, dupl, good, pass1
    """
    th = evas_thresholds(extra_args)
    lo, hi = th["nf_borderline_lower"], th["nf_borderline_upper"]
    min_ew = th["no_strandedness_min_extreme_wf"]
    rec = pd.Series("", index=detail.index, dtype=object)
    why = pd.Series("", index=detail.index, dtype=object)
    if "evas_status" not in detail.columns:
        return rec, why

    d = detail.copy()
    if info is not None and len(info):
        inf = info[["cell", "dupl", "mapped", "pass1"]].copy()
        inf["mc_dupl_rate"] = inf["dupl"] / inf["mapped"].where(inf["mapped"] > 0)
        d = d.merge(inf.rename(columns={"pass1": "mc_pass1"}), on="cell", how="left")
        d.index = detail.index
    passing = d["evas_status"] == "pass"
    ipd_fences = (
        _fences(d.loc[passing, "evas_mean_interpoint_distance"])
        if "evas_mean_interpoint_distance" in d
        else None
    )
    dup_ref = (
        d.loc[passing, "mc_dupl_rate"].median()
        if "mc_dupl_rate" in d and passing.any()
        else np.nan
    )
    has_ashleys = "ashleys_probability" in d.columns
    keep_prob = max(0.8, float(ashleys_threshold))

    for i, r in d.loc[d["evas_status"] == "borderline"].iterrows():
        bad, weak, notes = [], [], []

        nf = r.get("evas_noisy_fraction")
        if pd.notna(nf) and hi > lo:
            pos = (nf - lo) / (hi - lo)
            third = "lower" if pos <= 1 / 3 else ("upper" if pos >= 2 / 3 else "middle")
            notes.append(f"noisy fraction {nf:.3f} ({third} third of {lo:g}-{hi:g})")
            if third != "lower":
                weak.append("noisy fraction")
        else:
            weak.append("noisy fraction")
            notes.append("noisy fraction missing")

        ew = r.get("evas_n_extreme_wf")
        if pd.notna(ew) and ew < min_ew:
            bad.append("strandedness")
            notes.append(
                f"n_extreme_wf {ew:.0f} < {min_ew:g} (would fail 'no strandedness')"
            )
        elif pd.notna(ew):
            notes.append(f"n_extreme_wf {ew:.0f} >= {min_ew:g}")
        else:
            weak.append("strandedness")

        tot = r.get("evas_tot_count")
        if pd.notna(tot) and tot >= 2 * th["min_reads"]:
            notes.append(f"{tot:,.0f} reads")
        else:
            weak.append("depth")
            notes.append(
                f"{tot:,.0f} reads (< {2 * th['min_reads']:,})"
                if pd.notna(tot)
                else "reads missing"
            )

        ipd = r.get("evas_mean_interpoint_distance")
        if ipd_fences is None or pd.isna(ipd):
            weak.append("interpoint distance")
            notes.append("interpoint distance not comparable (< 5 pass cells)")
        elif not ipd_fences[0] <= ipd <= ipd_fences[1]:
            weak.append("interpoint distance")
            notes.append("interpoint distance outside the pass cells' range")

        if has_ashleys:
            p = r.get("ashleys_probability")
            if pd.isna(p):
                weak.append("ashleys-qc")
            elif p < float(ashleys_threshold):
                bad.append("ashleys-qc")
                notes.append(f"ashleys-qc {p:.2f} < {ashleys_threshold:g}")
            elif p < keep_prob:
                weak.append("ashleys-qc")
                notes.append(f"ashleys-qc {p:.2f} (< {keep_prob:g})")
            else:
                notes.append(f"ashleys-qc {p:.2f}")
        else:
            weak.append("ashleys-qc")
            notes.append("ashleys-qc not run")

        if "mc_pass1" in d.columns and pd.notna(r.get("mc_pass1")):
            if r["mc_pass1"] == 0:
                bad.append("coverage")
                notes.append("mosaicatcher pass1 = 0")
            dr = r.get("mc_dupl_rate")
            if pd.notna(dr) and pd.notna(dup_ref) and dr > dup_ref + 0.10:
                weak.append("duplicates")
                notes.append(f"duplicate rate {dr:.2f} (pass cells {dup_ref:.2f})")
            elif pd.notna(dr):
                notes.append(f"duplicate rate {dr:.2f}")
        else:
            weak.append("MosaiCatcher counts")
            notes.append("MosaiCatcher counts missing")

        if bad:
            rec[i] = "discard"
            head = "bad: " + ", ".join(bad)
        elif weak:
            rec[i] = "review"
            head = "not all signals agree: " + ", ".join(weak)
        else:
            rec[i] = "keep"
            head = "all signals agree"
        why[i] = f"{head} | " + "; ".join(notes)
    return rec, why
