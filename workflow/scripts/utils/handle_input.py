import pandas as pd
import os, sys


class HandleInput:
    def __init__(self, input_path, output_path, check_sm_tag=False, bam=True):
        df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_config_files.to_csv(output_path, sep="\t", index=False)
        self.df_config_files = df_config_files

    @staticmethod
    def handle_input_data(thisdir, exclude_list=list, bam=bool):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        ext = ".bam" if bam is True else ".fastq.gz"
        folder = "all" if bam is True else "fastq"
        complete_df_list = list()
        # print(thisdir)
        for sample in [e for e in os.listdir(thisdir) if e not in ["config", "log", ".DS_Store", "._.DS_Store"]]:
            # print(thisdir, sample, folder, ext)

            # print("{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder))
            l_files_all = [
                f
                for f in os.listdir("{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder))
                if f.endswith(ext)
            ]
            df = pd.DataFrame([{"File": f} for f in l_files_all])
            df["File"] = df["File"].str.replace(ext, "", regex=True)
            df["Folder"] = thisdir
            df["Sample"] = sample
            df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
            df["Full_path"] = "{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder)
            df["Full_path"] = df["Full_path"] + df["File"] + ext
            # if bam is True:s
            # l_files_selected = [f for f in os.listdir(thisdir + "/" + sample + "/selected/") if f.endswith(".bam")]
            # print(l_files_selected)
            # join = list(set(l_files_all).intersection(set(l_files_selected)))
            # df["Selected"] = False
            # df.loc[df["File"].isin(join), "Selected"] = True

            complete_df_list.append(df)

        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(drop=True)
        # complete_df = complete_df.loc[~complete_df["Cell"].isin(exclude_list)]
        return complete_df
