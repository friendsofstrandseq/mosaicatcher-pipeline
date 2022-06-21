import pandas as pd
import sys, os

input_dir = sys.argv[1]
input_dir += "/" if input_dir.endswith("/") is False else input_dir

l_files_all = [f for f in os.listdir(input_dir + "all/") if f.endswith('.bam')]
l_files_selected = [f for f in os.listdir(input_dir + "selected/") if f.endswith('.bam')]
join = list(set(l_files_all).intersection(set(l_files_selected)))
df = pd.DataFrame([{"File" : f} for f in l_files_all])
df.loc[df["File"].isin(join), "Selected"] = True
df["Selected"] = df["Selected"].fillna(False)
