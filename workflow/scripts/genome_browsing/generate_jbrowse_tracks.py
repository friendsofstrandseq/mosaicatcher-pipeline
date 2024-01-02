import pyBigWig
import pandas as pd
import sys, os
import subprocess

# Assuming sv_cell_df already contains the 'sv_call_name' column


# Function to map sv_call_name to color
def map_color(sv_call_name):
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
    return colors.get(sv_call_name, "#000000")  # Default to black if not found


# Load read counts data
counts_file_init = pd.read_csv(sys.argv[1], sep="\t", compression="gzip")
chrom_size_df = counts_file_init.groupby("chrom")["end"].max().reset_index()

# Load SV data
sv_df = pd.read_csv(sys.argv[2], sep="\t")

output_file = sys.argv[3]
output_dir = "/".join(output_file.split("/")[:-1])

os.makedirs(output_dir, exist_ok=True)


# Process each cell
for cell in sorted(counts_file_init.cell.unique()):
    print(cell)
    # Filter read counts for current cell
    counts_file = counts_file_init[counts_file_init["cell"] == cell]
    print(counts_file)
    counts_file["w"] = counts_file["w"].astype(float)
    counts_file["c"] = counts_file["c"].astype(float) * -1

    # Create BigWig for Watson counts
    bw = pyBigWig.open(f"{output_dir}/{cell}-W.bigWig", "w")
    bw.addHeader(list(chrom_size_df.itertuples(index=False, name=None)))
    bw.addEntries(
        counts_file.chrom.values.tolist(),
        counts_file.start.values.tolist(),
        ends=counts_file.end.values.tolist(),
        values=counts_file.w.values.tolist(),
    )
    bw.close()

    # Create BigWig for Crick counts
    bw = pyBigWig.open(f"{output_dir}/{cell}-C.bigWig", "w")
    bw.addHeader(list(chrom_size_df.itertuples(index=False, name=None)))
    bw.addEntries(
        counts_file.chrom.values.tolist(),
        counts_file.start.values.tolist(),
        ends=counts_file.end.values.tolist(),
        values=counts_file.c.values.tolist(),
    )
    bw.close()

    # Process SV data for current cell
    sv_cell_df = sv_df[sv_df["cell"] == cell]
    # sv_cell_df = sv_cell_df[["chrom", "start", "end", "sv_call_name"]]

    # Add a color column
    sv_cell_df["color"] = sv_cell_df["sv_call_name"].apply(map_color)

    # Select the relevant columns for the BED file
    sv_cell_df = sv_cell_df[
        [
            "chrom",
            "start",
            "end",
            "sv_call_name",
            "color",
            "sv_call_haplotype",
            "llr_to_ref",
            "af",
        ]
    ]

    # Write to BED file
    sv_filename = f"{output_dir}/{cell}-SV.bed"
    sv_cell_df.to_csv(sv_filename, sep="\t", index=False, header=False)

    # Compress the file using bgzip
    compressed_filename = sv_filename + ".gz"
    subprocess.run(["bgzip", "-c", sv_filename], stdout=open(compressed_filename, "wb"))

    # Index the compressed file using tabix
    subprocess.run(["tabix", "-p", "bed", compressed_filename])

    # Create BigWig for SV - this step might need adjustment based on SV data format
    # bw = pyBigWig.open(sv_filename, "w")
    # bw.addHeader(list(chrom_size_df.itertuples(index=False, name=None)))
    # Additional logic for adding SV data to BigWig may be required here
    # bw.close()

from pathlib import Path

Path(output_file).touch()
