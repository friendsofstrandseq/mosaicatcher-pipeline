import math
from collections import defaultdict

configfile: "config/Snake.config_embl.yaml"
import pandas as pd
import os, sys
from pprint import pprint
import pysam
from tqdm import tqdm

# TODO I/O : Function to define inputs ; simplify list/dict system // SOLVED
# TODO Use remote file system to download example files


def handle_input_data(thisdir, exclude_list=list):
    """ """
    # Parsing folder and retrieve only files with .bam extension
    data = [(r, file.replace(".bam", "")) for r, d, f in os.walk(thisdir) for file in f if ".bam" in file and ".bai" not in file]

    # Building pandas df based on folder structure
    df = pd.DataFrame(data, columns=["Folder", "File"])

    # Defining cols
    df["all/selected"] = df["Folder"].apply(lambda r: r.split("/")[-1])
    df["Sample"] = df["Folder"].apply(lambda r: r.split("/")[-2])
    df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
    df["Full_path"] = df["Folder"] + "/" + df["File"] + ".bam"

    # Filtering based on exclude list defined
    df_config_files = df.loc[~df["Cell"].isin(exclude_list)]

    # Export dicts
    SAMPLES = sorted(df_config_files.Sample.unique().tolist())
    BAM_PER_SAMPLE = df_config_files.loc[df_config_files["all/selected"] == "selected"].groupby("Sample")["File"].apply(list).to_dict()
    CELL_PER_SAMPLE = df_config_files.loc[df_config_files["all/selected"] == "selected"].groupby("Sample")["Cell"].apply(list).to_dict()
    ALLBAMS_PER_SAMPLE = df_config_files.loc[df_config_files["all/selected"] == "all"].groupby("Sample")["File"].apply(list).to_dict()

    return SAMPLES, BAM_PER_SAMPLE, CELL_PER_SAMPLE, ALLBAMS_PER_SAMPLE, df_config_files


def check_bam_header(bam_file_path):
    """ """

    # Get BAM file header with pysam
    h = pysam.view("-H", bam_file_path)
    h = [e.split("\t") for e in h.split("\n")]
    sm_tag_list = list(set([sub_e.replace("SM:", "") for e in h for sub_e in e if "SM:" in sub_e]))

    # Folder name based on path
    folder_name = bam_file_path.split("/")[-3]

    # Assertions
    assert len(sm_tag_list) == 1, "Two different SM tags in the header of BAM file {}".format(bam_file_path)
    assert sm_tag_list[0] == folder_name, 'Folder name "{}" must correspond to SM tag in BAM file "{}"'.format(folder_name, bam_file_path)


# FIXME : tmp solution to remove bad cells => need to fix this with combination of ASHLEYS ?
# TODO : other solution by giving in config file, CLI input ?

exclude_list = ["BM510x3PE20490"]

SAMPLES, BAM_PER_SAMPLE, CELL_PER_SAMPLE, ALLBAMS_PER_SAMPLE, df_config_files = handle_input_data(
    thisdir="/g/korbel2/weber/MosaiCatcher_files/bam_KG_full/", exclude_list=exclude_list
)

# print(df_config_files)
# print(df_config_files["Full_path"][0])

tqdm.pandas(desc="Checking if BAM SM tags correspond to folder names")
df_config_files["Full_path"].progress_apply(
    check_bam_header,
)

all_dict = df_config_files.loc[df_config_files["all/selected"] == "all"].groupby("Sample")["Cell"].nunique().to_dict()
selected_dict = df_config_files.loc[df_config_files["all/selected"] == "selected"].groupby("Sample")["Cell"].nunique().to_dict()
print("Detected {} samples:".format(df_config_files.Sample.nunique()))
[print("  {}:\t{} cells\t {} selected cells".format(s, all_dict[s], selected_dict[s])) for s in SAMPLES]

