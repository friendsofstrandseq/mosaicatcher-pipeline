import pandas as pd
import pysam 
import os, sys
import subprocess
from io import StringIO
from tqdm import tqdm
import parmap
import multiprocessing as mp

dirs = [
    "/g/korbel2/weber/MosaiCatcher_files/PRJEB30027/bam/C7_data/selected/",
    "/g/korbel2/weber/MosaiCatcher_files/PRJEB30027/bam/RPE1-WT/selected/",
    "/g/korbel2/weber/MosaiCatcher_files/PRJEB30027/bam/RPE-BM510/selected/",
    "/g/korbel2/weber/MosaiCatcher_files/LCL/H2NCTAFX2_GM20509B_20s000579-1-1/selected/",
]

for directory_input in dirs: 
    l_files_selected = [f for f in os.listdir(directory_input) if f.endswith('.bam')]
    sample = directory_input.split("/")[-3]
    m = mp.Manager()
    l_df = m.list()

    # for file in tqdm(l_files_selected):
        # print(file)
    def loop(file, l_df):
        p = subprocess.Popen('samtools idxstats ' + directory_input + file, shell=True, stdout=subprocess.PIPE)
        # text = [sub_e.split(" ") for e in p.communicate()[0].decode('utf-8').strip().split('\n') for sub_e in e.lstrip().split('\t')]
        text = [e.split('\t') for e in p.communicate()[0].decode('utf-8').strip().split('\n')]

        df = pd.DataFrame.from_records(text, columns=["chr", "chr_length", "read-segments_mapped", "read-segments_unmapped"])
        df = df.loc[df["chr"].isin(["chrX", "chrY"])]
        df[["chr_length", "read-segments_mapped", "read-segments_unmapped"]] = df[["chr_length", "read-segments_mapped", "read-segments_unmapped"]].astype(int)
        df["Ratio"] = df["read-segments_mapped"] / df["chr_length"]
        xy_ratio = df.loc[df["chr"] == "chrX", "Ratio"].values[0] / df.loc[df["chr"] == "chrY", "Ratio"].values[0]
        new_df = pd.DataFrame([
                {
                    "chrX/chrY_ratio" : xy_ratio, 
                    "Cell" : file.replace(".sort.mdup.bam", ""),
                    "Sample" : sample
                }]
            )
        l_df.append(new_df)


    parmap.starmap(loop, list(zip(l_files_selected)), l_df, pm_pbar=True, pm_processes=10)

    final_df = pd.concat(list(l_df))
    final_df.loc[final_df["chrX/chrY_ratio"] >= 5, "M/F"] = "F"
    final_df["M/F"] = final_df["M/F"].fillna("M")

    print(final_df)

    mf_dict_from_df = final_df["M/F"].value_counts().to_dict()
    mf_dict_from_df.setdefault("M", 0)
    mf_dict_from_df.setdefault("F", 0)
    print(mf_dict_from_df)

    try:
        perc_m = (mf_dict_from_df["M"] - mf_dict_from_df["F"]) / mf_dict_from_df["M"]
        perc_f = 1 - perc_m
    except:
        perc_f = (mf_dict_from_df["F"] - mf_dict_from_df["M"]) / mf_dict_from_df["F"]
        perc_m = 1 - perc_f

    print(perc_f, perc_m)
    sex = "M" if perc_m > 0.5 else "F"
    print(sex)
    final_df.to_csv("test_chrxy_analysis{}.tsv".format(sample), sep="\t")

