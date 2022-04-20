import os, sys, shutil
import pandas as pd
from tqdm import tqdm
from pandarallel import pandarallel

pandarallel.initialize(nb_workers=30, progress_bar=True)

ena_dl_location = sys.argv[1]
# ena_dl_location = "/g/korbel2/weber/MosaiCatcher_files/PRJEB30027/"

ena_dl_location = ena_dl_location + "/" if ena_dl_location.endswith("/") is False else ena_dl_location
bam_folder_location = ena_dl_location + "bam/"

# ENA FOLDERS & SUBFOLDERS STRUCTURE
data = [(r, file, r + "/" + file) for r, d, f in os.walk(ena_dl_location) for file in f if "sort.mdup.bam" in file and "/bam/" not in r]

# BUILD PANDAS DF
df = pd.DataFrame(data, columns=["Folder", "File", "Full_path"])
df["Sample_raw"] = df["File"].apply(lambda r: r.split("PE2")[0].split("x")[0])

# REAL SM TAGS NAME CORRESPONDING TO SAMPLES
dict_map = {"RPE1WT": "RPE1-WT", "BM510": "RPE-BM510", "C7": "C7_data"}
df["Sample_SM_format"] = df["Sample_raw"].map(dict_map)

pd.options.display.max_colwidth = 150
print(df)

# CREATE DIRS USED BY MC
[os.makedirs(bam_folder_location + e + "/all/", exist_ok=True) for e in df["Sample_SM_format"]]
[os.makedirs(bam_folder_location + e + "/selected/", exist_ok=True) for e in df["Sample_SM_format"]]

## Single threaded
# tqdm.pandas(desc='Moving files to bam/"sample"/all/ folder')
# df.progress_apply(lambda r: shutil.copy(r["Full_path"], bam_folder_location + r["Sample_SM_format"] + "/all/" + r["File"]), axis=1)
# tqdm.pandas(desc='Moving files to bam/"sample"/selected/ folder')
# df.progress_apply(lambda r: shutil.copy(r["Full_path"], bam_folder_location + r["Sample_SM_format"] + "/selected/" + r["File"]), axis=1)

## Multi-threaded / also long
print('Moving files to bam/"sample"/all/ folder')
df.parallel_apply(lambda r: shutil.copy(r["Full_path"], bam_folder_location + r["Sample_SM_format"] + "/all/" + r["File"]), axis=1)
print('Moving files to bam/"sample"/selected/ folder')
df.parallel_apply(lambda r: shutil.copy(r["Full_path"], bam_folder_location + r["Sample_SM_format"] + "/selected/" + r["File"]), axis=1)
