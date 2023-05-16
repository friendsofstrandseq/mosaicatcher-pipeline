import pandas as pd
import os, sys, gzip

colors = {
    "none": "248,248,248",  # #F8F8F8
    "del_h1": "119,170,221",  # #77AADD
    "del_h2": "68,119,170",  # #4477AA
    "del_hom": "17,68,119",  # #114477
    "dup_h1": "204,153,187",  # #CC99BB
    "dup_h2": "170,68,136",  # #AA4488
    "dup_hom": "119,17,85",  # #771155
    "inv_h1": "221,221,119",  # #DDDD77
    "inv_h2": "170,170,68",  # #AAAA44
    "inv_hom": "119,119,17",  # #777711
    "idup_h1": "221,170,119",  # #DDAA77
    "idup_h2": "170,119,68",  # #AA7744
    "complex": "119,68,17",  # #774411
}


input_file = sys.argv[1]  # Replace with your input TSV file
output_file = sys.argv[2]  # Replace with your desired output BED file

# Read the input TSV file
df = pd.read_csv(input_file, sep="\t")
df["strand"] = "+"
df["score"] = 0

# Map sv_call_name to colors
df["color"] = df["sv_call_name"].map(colors)


# Sort the DataFrame alphabetically by the 'cell' column
df_bed = df.sort_values(by=["cell"])

# Save DataFrame to BED file with separate tracks for each cell
with gzip.open(output_file, "wt") as f:
    f.write("browser position chr17\n")

    # Iterate over unique cells and create a track for each one
    for cell in df_bed["cell"].unique():
        df_cell = df_bed[df_bed["cell"] == cell]

        # Reorder columns to comply with BED format (chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb)
        tmp_df_bed = df_cell[["chrom", "start", "end", "sv_call_name", "score", "strand", "start", "end", "color"]]

        f.write(f'track name="{cell}" description="Track for cell {cell}" visibility=squish itemRgb="On"\n')
        tmp_df_bed.to_csv(f, sep="\t", header=False, index=False, mode="a", compression="gzip")
