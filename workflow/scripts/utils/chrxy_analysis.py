import pandas as pd
import pysam 
import os, sys
import subprocess
from io import StringIO
from tqdm import tqdm
import parmap
import multiprocessing as mp

# Snakemake input 
directory_input = snakemake.input[0]
directory_input += "/" if directory_input.endswith("/") is False else directory_input

# List bam files
l_files_selected = [f for f in os.listdir(directory_input) if f.endswith('.bam')]

# Retrieve sample name
sample = directory_input.split("/")[-3]

# Initiate MP
m = mp.Manager()
l_df = m.list()


def loop(file, l_df):
    """MP function to compute chrX/chrY coverage ratio 

    Args:
        file (str): bam file path
        l_df (list): MP shared list 
    """
    # Samtools shell command script
    p = subprocess.Popen('samtools idxstats ' + directory_input + file, shell=True, stdout=subprocess.PIPE)

    # Retrieve text from stdout & process it in python
    text = [e.split('\t') for e in p.communicate()[0].decode('utf-8').strip().split('\n')]

    # Convert list of list to pandas dataframe
    df = pd.DataFrame.from_records(text, columns=["chr", "chr_length", "read-segments_mapped", "read-segments_unmapped"])

    # Filter to keep chrX & Y
    df = df.loc[df["chr"].isin(["chrX", "chrY"])]

    # Convert as int type
    df[["chr_length", "read-segments_mapped", "read-segments_unmapped"]] = df[["chr_length", "read-segments_mapped", "read-segments_unmapped"]].astype(int)

    # Compute ratio
    df["Ratio"] = df["read-segments_mapped"] / df["chr_length"]
    xy_ratio = df.loc[df["chr"] == "chrX", "Ratio"].values[0] / df.loc[df["chr"] == "chrY", "Ratio"].values[0]

    # Create new df
    new_df = pd.DataFrame([
            {
                "chrX/chrY_ratio" : xy_ratio, 
                "Cell" : file.replace(".sort.mdup.bam", ""),
                "Sample" : sample
            }]
        )
    # Add to shared MP list
    l_df.append(new_df)

# Launch function in parallel on list of files
parmap.starmap(loop, list(zip(l_files_selected)), l_df, pm_pbar=True, pm_processes=10)

# Concatenate & Identify Male & Female cells based on cutoff
cutoff = 5
final_df = pd.concat(list(l_df))
final_df.loc[final_df["chrX/chrY_ratio"] >= cutoff, "M/F"] = "F"
final_df["M/F"] = final_df["M/F"].fillna("M")

print(final_df)

# Snakemake output
final_df.to_csv(snakemake.output.sex_analysis_cellwise, sep="\t", index=False)

mf_dict_from_df = final_df["M/F"].value_counts().to_dict()
mf_dict_from_df.setdefault("M", 0)
mf_dict_from_df.setdefault("F", 0)

try:
    perc_m = (mf_dict_from_df["M"] - mf_dict_from_df["F"]) / mf_dict_from_df["M"]
    perc_f = 1 - perc_m
except:
    perc_f = (mf_dict_from_df["F"] - mf_dict_from_df["M"]) / mf_dict_from_df["F"]
    perc_m = 1 - perc_f

sex = "M" if perc_m > 0.5 else "F"

with open(snakemake.output.sex_analysis_samplewise, "w") as output:
    output.write("{}\t{}".format(sample, sex))
