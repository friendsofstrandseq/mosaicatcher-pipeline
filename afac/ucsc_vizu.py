import pandas as pd
import os, sys, glob, gzip

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


def create_bed_row(row, category, color):
    chrom, start, end, cell = row["chrom"], row["start"], row["end"], row["cell"]
    score = "0"
    strand = "+"
    return f"{chrom}\t{start}\t{end}\t{category}\t{score}\t{strand}\t{start}\t{end}\t{color}\n"


def process_file(input_file, df_sv, output):
    # Extract cell name
    cell_name = os.path.basename(input_file).replace(".txt.percell.gz", "")

    # Read the input gzipped file
    df = pd.read_csv(input_file, sep="\t")

    # Create separate DataFrames for 'c' and 'w' columns
    df_c = df[["chrom", "start", "end", "c"]]
    df_c["c"] = df_c["c"] * -1
    df_w = df[["chrom", "start", "end", "w"]]

    # Filter df_sv
    df_sv_cell = df_sv.loc[df_sv["cell"] == cell_name]

    with gzip.open(output, "at") as output_file:
        output_file.write(
            f"track type=bedGraph name={cell_name}_W maxHeightPixels=40 description=BedGraph_{cell_name}_w.sort.mdup.bam_allChr visibility=full color=244,163,97\n"
        )
        df_w.to_csv(output_file, compression="gzip", sep="\t", header=False, index=False, mode="a")

        output_file.write(
            f"track type=bedGraph name={cell_name}_C maxHeightPixels=40 description=BedGraph_{cell_name}_c.sort.mdup.bam_allChr visibility=full color=102,139,138\n"
        )
        df_c.to_csv(output_file, compression="gzip", sep="\t", header=False, index=False, mode="a")

        output_file.write(f'track name="{cell_name}_SV" description="SV_call_name for cell {cell_name}" visibility=squish itemRgb="On"\n')
        for _, row in df_sv_cell.iterrows():
            bed_row = create_bed_row(row, row["sv_call_name"], row["color"])
            output_file.write(bed_row)


def main(input_counts_folder, input_sv_file, output_file):
    # Concatenate the W, C, and SV_call_name DataFrames
    df_sv = pd.read_csv(input_sv_file, sep="\t")  # Replace with your input TSV file containing SV_call_name information
    df_sv["color"] = df_sv["sv_call_name"].map(colors)
    df_sv = df_sv.sort_values(by=["cell"])

    # Get the list of input files in the input folder
    input_files = glob.glob(os.path.join(input_counts_folder, "*.txt.percell.gz"))

    # Process each input file
    for input_file in sorted(input_files):
        process_file(input_file, df_sv, output_file)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_counts_folder>  <input_sv_file> <output_file>")
        sys.exit(1)

    input_counts_folder = sys.argv[1]
    input_sv_file = sys.argv[2]
    output_file = sys.argv[3]
    main(input_counts_folder, input_sv_file, output_file)
