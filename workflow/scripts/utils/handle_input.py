import pandas as pd
import os, sys

# Simple class to retrieve automatically files in the fastq/bam folder and create a config dataframe
class HandleInput:
    def __init__(
        self,
        input_path,
        output_path,
        check_sm_tag=False,
        bam=True,
        genecore=False,
        genecore_path=str,
    ):
        if genecore is False:
            df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
        elif genecore is True:
            df_config_files, d_master = self.handle_input_data_genecore(thisdir=genecore_path)
            self.d_master = d_master

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_config_files.to_csv(output_path, sep="\t", index=False)
        self.df_config_files = df_config_files

    @staticmethod
    def handle_input_data_genecore(thisdir):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        complete_df_list = list()

        # List of folders/files to not consider (restrict to samples only)
        l = [
            e
            for e in os.listdir(
                "{genecore_prefix}/{date_folder}".format(
                    genecore_prefix=config["genecore_prefix"],
                    date_folder=config["genecore_date_folder"],
                )
            )
            if e.endswith(".gz")
        ]

        # Create a list of  files to process for each sample
        d_master = collections.defaultdict(dict)
        sub_l = list()
        for j, e in enumerate(l):
            sub_l.append(e)
            if (j + 1) % 192 == 0:
                common_element = findstem(sub_l)
                l_elems = common_element.split("lane1")
                prefix = l_elems[0]
                technician_name = l_elems[0].split("_")[-2]
                sample = l_elems[1].split("x")[0]
                index = l_elems[1].split("x")[1].split("PE")[0][-1]
                pe_index = common_element[-1]
                sub_l = list()

                d_master[sample]["prefix"] = prefix
                d_master[sample]["technician_name"] = technician_name
                d_master[sample]["index"] = index
                d_master[sample]["pe_index"] = pe_index
                d_master[sample]["common_element"] = common_element

        samples_to_process = config["samples_to_process"] if len(config["samples_to_process"]) > 0 else list(d_master.keys())

        config["data_location"] = "{data_location}/{genecore_date_folder}".format(
            data_location=config["data_location"],
            genecore_date_folder=config["genecore_date_folder"],
        )

        genecore_list = [
            expand(
                "{data_location}/{sample}/fastq/{sample}x0{index}PE20{cell_nb}.{pair}.fastq.gz",
                data_location=config["data_location"],
                sample=sample,
                index=d_master[sample]["index"],
                cell_nb=list(
                    range(
                        (int(d_master[sample]["pe_index"]) * 100) + 1,
                        (int(d_master[sample]["pe_index"]) * 100) + 97,
                    )
                ),
                pair=["1", "2"],
            )
            for sample in d_master
            if sample in samples_to_process
        ]
        genecore_list = [sub_e for e in genecore_list for sub_e in e]

        complete_df_list = list()

        for sample in d_master:
            df = pd.DataFrame([{"File": os.path.basename(f), "Folder": os.path.dirname(f)} for f in genecore_list if sample in f])
            df["File"] = df["File"].str.replace(".fastq.gz", "", regex=True)
            df["Sample"] = sample
            df["Pair"] = df["File"].apply(lambda r: r.split(".")[1])
            df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
            df["Full_path"] = df["Folder"] + "/" + df["File"] + ".fastq.gz"
            df["Genecore_path"] = (
                config["genecore_prefix"]
                + "/"
                + config["genecore_date_folder"]
                + "/"
                + d_master[sample]["prefix"]
                + "lane1"
                + df["File"].str.replace(".", "_", regex=True)
                + "_sequence.txt.gz"
            )
            df["Genecore_file"] = d_master[sample]["prefix"] + "lane1" + df["File"].str.replace(".", "_", regex=True)
            df["Genecore_file"] = df["Genecore_file"].apply(lambda r: "_".join(r.split("_")[:-1]))

            # Concat dataframes for each sample & output
            complete_df_list.append(df)

        complete_df = pd.concat(complete_df_list)

        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(drop=True)
        pd.options.display.max_colwidth = 200
        print(complete_df)
        # exit()
        return complete_df, d_master

    @staticmethod
    def handle_input_data(thisdir, exclude_list=list, bam=bool):
        """_summary_

        Args:
            thisdir (_type_): _description_
            exclude_list (_type_, optional): _description_. Defaults to list.

        Returns:
            _type_: _description_
        """
        # Extension & folder based on bam boolean input
        ext = ".bam" if bam is True else ".fastq.gz"
        folder = "bam" if bam is True else "fastq"
        complete_df_list = list()
        # List of folders/files to not consider (restrict to samples only)
        exclude = [
            "._.DS_Store",
            ".DS_Store",
            "all",
            "ashleys_counts",
            "bam",
            "cell_selection",
            "config",
            "counts",
            "fastq",
            "fastqc",
            "haplotag",
            "log",
            "merged_bam",
            "mosaiclassifier",
            "normalizations",
            "ploidy",
            "plots",
            "predictions",
            "segmentation",
            "snv_calls",
            "stats",
            "strandphaser",
        ]

        for sample in [e for e in os.listdir(thisdir) if e not in exclude]:
            # Create a list of  files to process for each sample
            l_files_all = [
                f
                for f in os.listdir("{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder))
                if f.endswith(ext)
            ]

            # Dataframe creation
            df = pd.DataFrame([{"File": f} for f in l_files_all])
            df["File"] = df["File"].str.replace(ext, "", regex=True)
            df["Folder"] = thisdir
            df["Sample"] = sample
            df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
            df["Full_path"] = "{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder)
            df["Full_path"] = df["Full_path"] + df["File"] + ext

            complete_df_list.append(df)

        # Concat dataframes for each sample & output
        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(drop=True)
        return complete_df


# GENECORE


def findstem(arr):

    # Determine size of the array
    n = len(arr)

    # Take first word from array
    # as reference
    s = arr[0]
    l = len(s)

    res = ""

    for i in range(l):
        for j in range(i + 1, l + 1):

            # generating all possible substrings
            # of our reference string arr[0] i.e s
            stem = s[i:j]
            k = 1
            for k in range(1, n):

                # Check if the generated stem is
                # common to all words
                if stem not in arr[k]:
                    break

            # If current substring is present in
            # all strings and its length is greater
            # than current result
            if k + 1 == n and len(res) < len(stem):
                res = stem

    return res


# class HandleInput:
#     def __init__(self, input_path, output_path, check_sm_tag=False, bam=True):
#         df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
#         os.makedirs(os.path.dirname(output_path), exist_ok=True)
#         df_config_files.to_csv(output_path, sep="\t", index=False)
#         self.df_config_files = df_config_files

#     @staticmethod
#     def handle_input_data(thisdir, exclude_list=list, bam=bool):
#         """_summary_

#         Args:
#             thisdir (_type_): _description_
#             exclude_list (_type_, optional): _description_. Defaults to list.

#         Returns:
#             _type_: _description_
#         """
#         ext = ".bam" if bam is True else ".fastq.gz"
#         folder = "bam" if bam is True else "fastq"
#         complete_df_list = list()
#         exclude = [
#             "._.DS_Store",
#             ".DS_Store",
#             "all",
#             "ashleys_counts",
#             "bam",
#             "cell_selection",
#             "config",
#             "counts",
#             "fastq",
#             "fastqc",
#             "haplotag",
#             "log",
#             "merged_bam",
#             "mosaiclassifier",
#             "normalizations",
#             "ploidy",
#             "plots",
#             "predictions",
#             "segmentation",
#             "snv_calls",
#             "stats",
#             "strandphaser",
#         ]
#         for sample in [e for e in os.listdir(thisdir) if e not in exclude]:
#             l_files_all = [
#                 f
#                 for f in os.listdir("{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder))
#                 if f.endswith(ext)
#             ]
#             df = pd.DataFrame([{"File": f} for f in l_files_all])
#             df["File"] = df["File"].str.replace(ext, "", regex=True)
#             df["Folder"] = thisdir
#             df["Sample"] = sample
#             df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
#             df["Full_path"] = "{thisdir}/{sample}/{folder}/".format(thisdir=thisdir, sample=sample, folder=folder)
#             df["Full_path"] = df["Full_path"] + df["File"] + ext

#             complete_df_list.append(df)

#         complete_df = pd.concat(complete_df_list)
#         complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(drop=True)
#         return complete_df
