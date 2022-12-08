import pandas as pd
import pysam
import os, sys
from tqdm import tqdm

def check_bam_header(bam_file_path):
    """ """

    # Get BAM file header with pysam
    h = pysam.view("-H", bam_file_path)
    h = [e.split("\t") for e in h.split("\n") if e.startswith("@RG")]
    sm_tag = list(set([sub_e.replace("SM:", "") for e in h for sub_e in e if sub_e.startswith("SM:")]))
    id_tag = list(set([sub_e.replace("ID:", "") for e in h for sub_e in e if sub_e.startswith("ID:")]))

    # Folder name based on path
    folder_name = bam_file_path.split("/")[-3]
    file_name = bam_file_path.split("/")[-1].replace(".sort.mdup.bam", "")
    try:
        # Assertions
        assert len(sm_tag) == 1, "Two different SM tags in the header of BAM file {}".format(bam_file_path)
        assert len(id_tag) == 1, "Two different ID tags in the header of BAM file {}".format(bam_file_path)
        assert sm_tag[0] == folder_name, 'Folder name "{}" must correspond to SM tag in BAM file "{}"'.format(folder_name, bam_file_path)
        assert id_tag[0] == file_name, 'File name "{}" must correspond to SM tag in BAM file "{}"'.format(file_name, bam_file_path)
        return True
    except AssertionError:
        return False
        
        


# df_config_files = pd.read_csv(".tests/data_CHR17_NEW_BAM/config/config_df.tsv", sep="\t")
# print(df_config_files.Full_path.tolist())
l = ['.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20301.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20302.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20303.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20304.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20305.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20306.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20307.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20308.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20309.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20310.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20311.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20312.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20313.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20314.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20316.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20317.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20318.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20319.sort.mdup.bam', '.tests/data_CHR17_NEW_BAM/RPE-BM510/bam/BM510x04_PE20320.sort.mdup.bam']
l = snakemake.input.bam
df = pd.DataFrame(pd.Series(l, name="Full_path"))
df["Sample"] = df["Full_path"].apply(lambda r: r.split("/")[-3])
df["Cell"] = df["Full_path"].apply(lambda r: r.split("/")[-1].replace(".sort.mdup.bam", ""))
print(df)
exit()

tqdm.pandas(desc="Checking if BAM SM tags correspond to folder names")
df_config_files["assert"] = df_config_files["Full_path"].progress_apply(check_bam_header)
print(df_config_files)

snakemake.to_csv(output[0], sep="\t", index=False)

# os.makedirs(os.path.dirname(output_path), exist_ok=True)