#!/usr/bin/env python3
from ipysheet import sheet, column, to_dataframe
import ipywidgets as w
import pandas as pd
import numpy as np


def get_ipysheet(cell_per_sample, sample, mosaic_info, ashleys_labels):

    # Retrive cell list from dictionnary
    cell_list = cell_per_sample[sample]
    row_nb = len(cell_list)

    # Draw table
    s = sheet(
        rows=row_nb, columns=7, column_headers=["cell", "selected?", "mosaic pass", "ashleys pred", "ashleys proba", "Good reads", "%dupl"]
    )
    s.column_width = [6, 2, 2, 2, 2, 3, 2]
    s.layout = w.Layout(width="800px", height="100%")

    # Retrieve low-coverage information from mosaic count "info" file
    mosaic_info_df = pd.read_csv(mosaic_info, sep="\t", skiprows=13)
    print(mosaic_info_df)

    # mosaic_info_df.loc[~np.isfinite(mosaic_info_df["dupl"])]["dupl"]

    mosaic_info_df["%dupl"] = 100 * (
        mosaic_info_df.loc[np.isfinite(mosaic_info_df["dupl"])]["dupl"]
        / mosaic_info_df.loc[np.isfinite(mosaic_info_df["mapped"])]["mapped"]
    )
    mosaic_info_df["%dupl"] = mosaic_info_df["%dupl"].fillna(0)
    # mosaic_info_df["%dupl"] = 100 * (mosaic_info_df["dupl"] / mosaic_info_df["mapped"])
    mosaic_info_df["%dupl"] = mosaic_info_df["%dupl"].round(0)
    mosaic_info_df["good"] = mosaic_info_df["good"].astype(str)
    print(mosaic_info_df)

    mosaic_info_df["%dupl"] = mosaic_info_df["%dupl"].astype(int).astype(str)
    mosaic_info_df["pass1"] = mosaic_info_df["pass1"].astype(bool)
    mosaic_info_cells = mosaic_info_df["pass1"].values.tolist()
    # print(mosaic_info_df)

    # ashleys labels
    ashleys_labels = pd.read_csv(ashleys_labels, sep="\t").sort_values(by="cell")
    ashleys_labels["cell"] = ashleys_labels["cell"].str.replace(".sort.mdup.bam", "")

    # Final pred
    final_df = pd.merge(mosaic_info_df, ashleys_labels, on="cell")
    final_df.loc[(final_df["prediction"] == 1) & (final_df["pass1"] == True), "final_pred"] = 1
    final_df["final_pred"] = final_df["final_pred"].fillna(0)
    final_df["final_pred"] = final_df["final_pred"].astype(bool)
    final_pred = final_df["final_pred"].values.tolist()

    # Cleaning
    ashleys_labels["prediction"] = ashleys_labels["prediction"].replace({1: "✅", 0: "❌"})
    ashleys_predictions = ashleys_labels["prediction"].values.tolist()
    ashleys_labels["probability"] = ashleys_labels["probability"].round(3)
    ashleys_proba = ashleys_labels["probability"].values.tolist()

    # Fill columns
    column(0, cell_list)
    column(1, mosaic_info_cells)
    column(2, mosaic_info_df["pass1"].replace({True: "✅", False: "❌"}).values.tolist())
    column(3, ashleys_predictions)
    column(4, ashleys_proba)
    column(5, mosaic_info_df["good"].values.tolist())
    column(6, mosaic_info_df["%dupl"].values.tolist())
    return s
