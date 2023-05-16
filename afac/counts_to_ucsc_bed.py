import pandas as pd
import os, sys, glob, gzip


def process_file(input_file, output):
    # Extract cell name
    cell_name = os.path.basename(input_file).replace(".txt.percell.gz", "")

    # Read the input gzipped file
    df = pd.read_csv(input_file, sep="\t")

    # Create separate DataFrames for 'c' and 'w' columns
    df_c = df[["chrom", "start", "end", "c"]]
    df_c["c"] = df_c["c"] * -1
    df_w = df[["chrom", "start", "end", "w"]]

    # Save 'c' DataFrame to the output file
    additional_row_c = f"track type=bedGraph name={cell_name} maxHeightPixels=40 description=BedGraph_{cell_name}_c.sort.mdup.bam_allChr visibility=full color=102,139,138\n"
    additional_row_w = f"track type=bedGraph name={cell_name} maxHeightPixels=40 description=BedGraph_{cell_name}_w.sort.mdup.bam_allChr visibility=full color=244,163,97\n"

    with gzip.open(output, "at") as output_file:
        output_file.write(additional_row_w)
        df_w.to_csv(output_file, compression="gzip", sep="\t", header=False, index=False, mode="a")

        output_file.write(additional_row_c)
        df_c.to_csv(output_file, compression="gzip", sep="\t", header=False, index=False, mode="a")


def main(input_folder, outputfile):
    # Get the list of input files in the input folder
    input_files = glob.glob(os.path.join(input_folder, "*.txt.percell.gz"))

    # Process each input file
    for input_file in sorted(input_files):
        process_file(input_file, outputfile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_file>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_file = sys.argv[2]
    main(input_folder, output_file)
