import os, sys
import subprocess
from tqdm import tqdm
import parmap
import multiprocessing as mp
import numpy as np


l_files_selected = snakemake.input.bam


# Initiate MP
m = mp.Manager()
l_df = m.list()


def loop(file, l_df):
    """MP function to compute chrX/chrY coverage ratio

    Args:
        file (str): bam file path
        l_df (list): MP shared list
    """
    p = subprocess.Popen("samtools view -c -f 1 {file}".format(file=file), shell=True, stdout=subprocess.PIPE)
    p_out = int(p.communicate()[0].decode("utf-8").strip())

    # Add to shared MP list
    l_df.append(p_out)


# Launch function in parallel on list of files
parmap.starmap(loop, list(zip(l_files_selected)), l_df, pm_pbar=True, pm_processes=10)

l_df = list(l_df)

paired_end = True
if np.mean(l_df) == 0:
    paired_end = False
elif np.mean(l_df) == 0:
    if 0 in l_df:
        sys.exit("Mix of single-end and paired-end files")
    else:
        paired_end = True

print("Paired-end: {paired_end} for sample: {sample}".format(paired_end=paired_end, sample=snakemake.wildcards.sample))
print(l_files_selected)
print(l_df)

with open(snakemake.output.single_paired_end_detect, "w") as output:
    output.write(str(paired_end).upper())
