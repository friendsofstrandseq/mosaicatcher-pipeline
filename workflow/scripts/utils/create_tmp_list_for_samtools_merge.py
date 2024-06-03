# Open the output file in append mode
with open(snakemake.output[0], "a") as f:
    for cell in snakemake.input.bam:
        # Write each filename followed by a newline to the output file
        f.write(cell + "\n")
