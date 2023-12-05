import pandas as pd
import os


def main(labels_file, subclone_file, selected_folder, output_file):
    # Read labels.tsv
    labels_df = pd.read_csv(labels_file, sep="\t")
    labels_cells = set(
        labels_df["cell"].str.replace(".sort.mdup.bam", "").values.tolist()
    )

    # Read input_subclonality.txt
    input_subclonality = pd.read_csv(subclone_file, sep="\t")
    subclone_cells = set(input_subclonality["Filename"].values.tolist())

    # List files in selected/ folder and process filenames
    selected_cells = set(
        file.replace(".sort.mdup.bam", "")
        for file in os.listdir(selected_folder)
        if file.endswith(".sort.mdup.bam")
    )

    # Compare sets
    if labels_cells == subclone_cells == selected_cells:
        result = "PASS: All cell lists match."
    else:
        result = "FAIL: Cell lists do not match."

    # Logging details of the mismatch
    with open(output_file, "w") as output:
        output.write("Labels cells: {}\n".format(labels_cells))
        output.write("Subclone cells: {}\n".format(subclone_cells))
        output.write("Selected cells: {}\n".format(selected_cells))
        output.write("Discrepancy details:\n")
        output.write(
            "In labels but not in subclone: {}\n".format(labels_cells - subclone_cells)
        )
        output.write(
            "In subclone but not in labels: {}\n".format(subclone_cells - labels_cells)
        )
        output.write(
            "In labels but not in selected: {}\n".format(labels_cells - selected_cells)
        )
        output.write(
            "In selected but not in labels: {}\n".format(selected_cells - labels_cells)
        )
        output.write(result)


if __name__ == "__main__":
    # Extracting Snakemake input variables
    labels_file = snakemake.input.labels
    subclone_file = snakemake.input.subclone_list
    selected_folder = snakemake.input.selected_cells
    output_file = snakemake.output[0]

    main(labels_file, subclone_file, selected_folder, output_file)
