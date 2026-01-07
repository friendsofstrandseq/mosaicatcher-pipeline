import pandas as pd
import pysam
import os, sys


# df_config_files = pd.read_csv(, sep="\t")

# READ BAM FILE HEADER OF FIRST BAM IN THE PANDAS DF
h = pysam.view("-H", snakemake.input.bam[0])
# h = pysam.view("-H", os.listdir(snakemake.input.bam + "selected")[0])
h = [e.split("\t") for e in h.split("\n") if "@SQ" in e]

# CONVERT TO PANDAS DF
df_h = pd.DataFrame(h, columns=["TAG", "Contig", "LN"])

# PROCESS CONTIGS
output_h = pd.DataFrame(df_h["Contig"].str.replace("SN:", ""))
print(snakemake.params["chroms"], type(snakemake.params["chroms"]))
output_h = output_h.loc[~output_h["Contig"].isin(snakemake.params["chroms"])]

# print(output_h)

# EXPORT
output_h["Contig"].to_csv(snakemake.output[0], index=False, sep="\t", header=False)
