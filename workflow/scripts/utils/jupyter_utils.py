#!/usr/bin/env python3
from ipysheet import sheet, column, to_dataframe
import ipywidgets as w
import pandas as pd


def get_ipysheet(cell_per_sample, sample, mosaic_info):

    # Retrive cell list from dictionnary
    cell_list = cell_per_sample[sample]
    row_nb = len(cell_list)

    # Draw table
    s = sheet(rows=row_nb, columns=2, column_headers=["cell", "selected?"])
    s.column_width = [8, 2]
    s.layout = w.Layout(width="400px", height="100%")

    # Retrieve low-coverage information from mosaic count "info" file
    mosaic_info_df = pd.read_csv(mosaic_info, sep="\t", skiprows=13)
    mosaic_info_df["pass1"] = mosaic_info_df["pass1"].astype(bool)
    mosaic_info_cells = mosaic_info_df["pass1"].values.tolist()

    # Fill columns
    column(0, cell_list)
    column(1, mosaic_info_cells)
    return s
