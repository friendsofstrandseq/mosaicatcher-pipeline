import sys, os
import pandas as pd

input_folder = snakemake.wildcards.input_folder
sample = snakemake.wildcards.sample

# input_folder = "/g/korbel2/weber/MosaiCatcher_files/HGSVC_WH"
# sample = "TEST"
ext = ".sort.mdup.bam"

# ASSERTIONS TO CHECK IF FOLDERS EXIST OR NOT
assert os.path.isdir(
    "{input_folder}/{sample}/bam/".format(input_folder=input_folder, sample=sample)
), "Folder all for sample {sample} does not exist".format(sample=sample)
assert os.path.isdir(
    "{input_folder}/{sample}/selected/".format(input_folder=input_folder, sample=sample)
), "Folder selected for sample {sample}  does not exist".format(sample=sample)

# RETRIEVE LIST OF FILES
l_files_all = [f for f in os.listdir("{input_folder}/{sample}/bam/".format(input_folder=input_folder, sample=sample)) if f.endswith(ext)]
l_files_selected = [
    f for f in os.listdir("{input_folder}/{sample}/selected/".format(input_folder=input_folder, sample=sample)) if f.endswith(ext)
]

# CHECK IF FILE EXTENSION IS CORRECT
if (
    len(l_files_all) == 0
    and len([f for f in os.listdir("{input_folder}/{sample}/bam/".format(input_folder=input_folder, sample=sample)) if f.endswith(".bam")])
    > 0
):
    sys.exit("BAM files extension were correctly set: .bam instead of .sort.mdup.bam (prevent further issues)")


# pd.options.display.max_rows = 100

# CREATE PANDAS DATAFRAME
df = pd.DataFrame([l_files_all, [1] * len(l_files_all), [1] * len(l_files_all)]).T
df.columns = ["cell", "probability", "prediction"]

# FLAG UNSELECTED FILES
unselected = set(l_files_all).difference(set(l_files_selected))
df.loc[df["cell"].isin(unselected), "probability"] = 0
df.loc[df["cell"].isin(unselected), "prediction"] = 0

# OUTPUT
df.to_csv(snakemake.output[0], sep="\t", index=False)
