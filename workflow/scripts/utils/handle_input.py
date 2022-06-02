from collections import defaultdict
import pandas as pd
import os, sys
from pprint import pprint
import pysam
from tqdm import tqdm


class HandleInput:
    def __init__(self, input_path, output_path, check_sm_tag):
        df_config_files = self.handle_input_data(thisdir=input_path)
      
        # if snakemake.config["check_sm_tag"] is True:
        # print(config["check_sm_tag"])
        # if config["check_sm_tag"] is True:
        if check_sm_tag is True:
            tqdm.pandas(desc="Checking if BAM SM tags correspond to folder names")
            df_config_files["Full_path"].progress_apply(self.check_bam_header)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_config_files.to_csv(output_path, sep="\t", index=False)
        self.df_config_files = df_config_files

    @staticmethod
    def handle_input_data(thisdir, exclude_list=list):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """

        complete_df_list = list()
        # input_dir = snakemake.config["input_bam_location"]
        # input_dir = config["input_bam_location"]
        # input_dir = input_path
        for sample in os.listdir(thisdir):
            l_files_all = [f for f in os.listdir(thisdir + sample + "/all/") if f.endswith('.bam')]
            l_files_selected = [f for f in os.listdir(thisdir + sample + "/selected/") if f.endswith('.bam')]
            join = list(set(l_files_all).intersection(set(l_files_selected)))
            df = pd.DataFrame([{"File" : f} for f in l_files_all])
            df.loc[df["File"].isin(join), "Selected"] = True
            df["File"] = df["File"].str.replace(".bam", "", regex=True)
            df["Selected"] = df["Selected"].fillna(False)
            df["Folder"] = thisdir
            df["Sample"] = sample
            df["Cell"] = df["File"].apply(lambda r : r.split(".")[0])
            df["Full_path"] = df["Folder"] + sample + "/all/" + df["File"] + ".bam"
            complete_df_list.append(df)
        complete_df = pd.concat(complete_df_list)

        # CHECKME : HGVSC errors
        exclude_list = ["GM18534Bx02PE20381", "HG02011x02PE20557", "HG02011x02PE20552"]
        
        
        complete_df = complete_df.loc[~complete_df["Cell"].isin(exclude_list)]
        return complete_df

    @staticmethod
    def check_bam_header(bam_file_path):
        """_summary_

        Args:
            bam_file_path (_type_): _description_
        """

        # Get BAM file header with pysam
        pysam.set_verbosity(0)
        h = pysam.view("-H", bam_file_path)
        h = [e.split("\t") for e in h.split("\n")]
        sm_tag_list = list(set([sub_e.replace("SM:", "") for e in h for sub_e in e if "SM:" in sub_e]))

        # Folder name based on path
        folder_name = bam_file_path.split("/")[-3]

        # Assertions
        assert len(sm_tag_list) == 1, "Two different SM tags in the header of BAM file {}".format(bam_file_path)
        assert sm_tag_list[0] == folder_name, 'Folder name "{}" must correspond to SM tag in BAM file "{}"'.format(
            folder_name, bam_file_path
        )

# if __name__ == "__main__":
#     c = HandleInput(snakemake.input[0], snakemake.output[0])
