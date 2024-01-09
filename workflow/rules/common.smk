import collections
import pandas as pd
import yaml

# from scripts.utils import handle_input, make_log_useful, pipeline_aesthetic_start
from scripts.utils import make_log_useful, pipeline_aesthetic_start
import os, sys

# Solve LC_CTYPE issue

os.environ["LC_CTYPE"] = "C"


# print(config["data_location"])

if config["ashleys_pipeline"] is True and config["genecore"] is True:
    config["data_location"] = config["abs_path"].join(
        config["data_location"].split("/")[:-1]
    )


envvars:
    "LC_CTYPE",


if config["list_commands"] is True:
    pipeline_aesthetic_start.argparse_help(config)


# Start with aesthetic pipeline config presentation
onstart:
    pipeline_aesthetic_start.pipeline_aesthetic_start(config)
    # subprocess.Popen(
    #     "rsync --ignore-existing -avzh config/config.yaml {folder_path}/config".format(
    #         folder_path=config["data_location"]
    #     ),
    #     shell=True,
    #     stdout=subprocess.PIPE,
    # )



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


if config["reference"] == "mm10":
    config["chromosomes"] = [
        "chr{e}".format(e=str(e)) for e in list(range(1, 20)) + ["X", "Y"]
    ]


if config["chromosomes_to_exclude"]:
    chroms_init = config["chromosomes"]
    chroms = [e for e in chroms_init if e not in config["chromosomes_to_exclude"]]
    config["chromosomes"] = chroms

# # List of assertions to verify
# if config["genecore"] is False:
#     l_samples = os.listdir(config["data_location"])
# #     assert (
# #         "fastq" not in l_samples
# #     ), "fastq folder found in the {} data_location specified: please specify a parent folder".format(
# #         config["data_location"]
# #     )

if config["ashleys_pipeline"] is True:
    assert (
        config["ashleys_pipeline"] != config["input_bam_legacy"]
    ), "ashleys_pipeline and input_bam_legacy parameters cannot both be set to True"


if (
    config["hgsvc_based_normalized_counts"] is True
    and config["multistep_normalisation_for_SV_calling"] is True
):
    assert (
        config["hgsvc_based_normalized_counts"]
        != config["multistep_normalisation_for_SV_calling"]
    ), "hgsvc_based_normalized_counts and multistep_normalisation_for_SV_calling parameters cannot both be set to True, parameters are mutually exclusive"


if config["multistep_normalisation_for_SV_calling"] is True:
    assert (
        config["multistep_normalisation_for_SV_calling"]
        == config["multistep_normalisation"]
    ), "multistep_normalisation parameter should be set to True"

# if (config["multistep_normalisation_for_SV_calling"] is True):
#     assert config["multistep_normalisation_for_SV_calling"] == config["ashleys_pipeline"], "ashleys_pipeline parameter should be set to True when multistep_normalisation_for_SV_calling is used"

if config["ashleys_pipeline_only"] is True:
    assert (
        config["ashleys_pipeline_only"] == config["ashleys_pipeline"]
    ), "ashleys_pipeline parameter should be set to True when ashleys_pipeline_only is also set to True"


if config["scNOVA"] is True:
    # print(config["chromosomes_to_exclude"])
    assert (
        "chrY" in config["chromosomes_to_exclude"]
    ), "chrY is not handled by scNOVA yet, please remove it for config['chromosomes'] and add it in config['chomosomes_to_exclude']"




if config["strandscape_labels_path"]:
    folder_location = config["abs_path"].join(
        config["strandscape_labels_path"].split("/")[:-1]
    )
    labels_path = f"{folder_location}/labels.tsv"
    assert os.path.isfile(labels_path)
    ashleys_labels = pd.read_csv(labels_path, sep="\t")
    strandscape_labels = pd.read_csv(config["strandscape_labels_path"], sep="\t")
    print(ashleys_labels)
    print(strandscape_labels)
    assert ashleys_labels.shape[0] == strandscape_labels.shape[0]
    assert (
        ashleys_labels.cell.values.tolist() == strandscape_labels.cell.values.tolist()
    )


# Configure if handle_input needs to be based on bam or fastq
bam = True if config["ashleys_pipeline"] is False else False


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
        # print(input_path)
        # print(genecore_path)
        # print("\n")
        if genecore is False:
            df_config_files = self.handle_input_data(thisdir=input_path, bam=bam)
        elif genecore is True:
            df_config_files, d_master = self.handle_input_data_genecore(
                thisdir=genecore_path
            )
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
        from pprint import pprint
        from collections import Counter

        directory_path = f"{config['genecore_prefix']}/{config['genecore_date_folder']}"

        l = sorted([e for e in os.listdir(directory_path) if e.endswith(".txt.gz")])

        complete_df_list = list()
        # print(thisdir)
        genecore_prefix = config["genecore_prefix"]
        date_folder = config["genecore_date_folder"]
        # print(f"{genecore_prefix}/{date_folder}")

        # Pattern to extract sample name and index
        pattern = re.compile(r"(.*_lane1)(.*?)(iTRU|PE20)(.*?)(\d{2})(?:_1_|_2_)")

        samples = list()
        prefixes = list()
        indexes = list()
        plate_types = list()
        d_master = collections.defaultdict(
            lambda: {
                "indexes": set(),
                "file_prefix": "",
                "plate_type": "",
                "index_pattern": "",
            }
        )

        # First pass: Count occurrences of each sample_name
        file_counts_per_sample = Counter()
        for file_path in l:
            match = pattern.search(file_path)
            if match:
                sample_name = match.group(2)
                file_counts_per_sample[sample_name] += 1

        # Second pass: Process files and determine plate type per sample
        for j, file_path in enumerate(sorted(l)):
            match = pattern.search(file_path)
            if match:
                sample_name = match.group(2)
                index = match.group(4)
                indexes.append(index)
                d_master[sample_name]["indexes"].add(index)
                file_count = file_counts_per_sample[sample_name]

                # Determine plate type using modulo 96 operation
                if file_count % 96 != 0:
                    raise ValueError(
                        f"Invalid file count for sample {sample_name} with file count {file_count}. Must be a multiple of 96."
                    )
                plate_type = int(file_count / 2)

                if (j + 1) % file_count == 0:
                    prefixes.append(match.group(3))
                    d_master[sample_name]["file_prefix"] = match.group(1)
                    d_master[sample_name]["index_pattern"] = match.group(3)
                    plate = directory_path.split("/")[-1]
                    samples.append(sample_name)
                    plate_types.append(plate_type)
                    d_master[sample_name]["plate_type"] = plate_type

        samples_to_process = (
            config["samples_to_process"]
            if len(config["samples_to_process"]) > 0
            else list(d_master.keys())
        )

        config["data_location"] = "{data_location}/{genecore_date_folder}".format(
            data_location=config["data_location"],
            genecore_date_folder=config["genecore_date_folder"],
        )

        genecore_list = [
            expand(
                "{data_location}/{sample}/fastq/{sample}{regex_element}{index}{cell_nb}.{pair}.fastq.gz",
                data_location=config["data_location"],
                sample=sample,
                regex_element=d_master[sample]["index_pattern"],
                index=d_master[sample]["indexes"],
                cell_nb=[str(e).zfill(2) for e in list(range(1, 97))],
                pair=["1", "2"],
            )
            for sample in d_master
            if sample in samples_to_process
        ]
        genecore_list = [sub_e for e in genecore_list for sub_e in e]
        # pprint(d_master)

        complete_df_list = list()

        for sample in d_master:
            df = pd.DataFrame(
                [
                    {"File": os.path.basename(f), "Folder": os.path.dirname(f)}
                    for f in genecore_list
                    if sample in f
                ]
            )
            if df.shape[0] > 0:
                df["File"] = df["File"].str.replace(".fastq.gz", "", regex=True)
                df["Sample"] = sample
                df["Pair"] = df["File"].apply(lambda r: r.split(".")[1])
                df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
                df["Full_path"] = df[["Folder", "File"]].apply(
                    lambda r: f"{r['Folder']}/{r['File']}.fastq.gz", axis=1
                )

                df["Genecore_path"] = df["File"].apply(
                    lambda r: f"{config['genecore_prefix']}/{config['genecore_date_folder']}/{d_master[sample]['file_prefix']}{r.replace('.', '_')}_sequence.txt.gz"
                )
                df["Genecore_file"] = df["File"].apply(
                    lambda r: f"{d_master[sample]['file_prefix']}{r.replace('.', '_')}"
                )
                df["Genecore_file"] = df["Genecore_file"].apply(
                    lambda r: "_".join(r.split("_")[:-1])
                )

                # Concat dataframes for each sample & output
                complete_df_list.append(df)

        complete_df = pd.concat(complete_df_list)

        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        pd.options.display.max_colwidth = 200
        # print(complete_df)
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

        l_to_process = [
            e
            for e in os.listdir(thisdir)
            if e not in exclude and e.endswith(".zip") is False
        ]
        # print(l_to_process)
        if config["samples_to_process"]:
            l_to_process = [
                e for e in l_to_process if e in config["samples_to_process"]
            ]

        for sample in l_to_process:
            # Create a list of  files to process for each sample
            l_files_all = [
                f
                for f in os.listdir(
                    "{thisdir}/{sample}/{folder}/".format(
                        thisdir=thisdir, sample=sample, folder=folder
                    )
                )
                if f.endswith(ext)
            ]
            # print(l_files_all)

            for f in l_files_all:
                if len(f.split("_")) == 4:
                    assert (
                        len(f.split("_")) != 4
                    ), "Your file name is using 4 times the '_' character, which is currently not supported by ashleys-qc, please rename your files"

            # Dataframe creation
            df = pd.DataFrame([{"File": f} for f in l_files_all])
            df["File"] = df["File"].str.replace(ext, "", regex=True)
            df["Folder"] = thisdir
            df["Sample"] = sample
            df["Cell"] = df["File"].apply(lambda r: r.split(".")[0])
            df["Full_path"] = "{thisdir}/{sample}/{folder}/".format(
                thisdir=thisdir, sample=sample, folder=folder
            )
            df["Full_path"] = df["Full_path"] + df["File"] + ext
            # print(df)

            complete_df_list.append(df)

        # Concat dataframes for each sample & output
        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
        # print(complete_df)
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


# Create configuration file with samples

# print("config['data_location']")
# print(config["data_location"])

c = HandleInput(
    input_path=config["data_location"],
    genecore_path="{genecore_prefix}".format(
        genecore_prefix=config["genecore_prefix"],
    ),
    # genecore_path="{genecore_prefix}/{genecore_date_folder}".format(
    #     genecore_prefix=config["genecore_prefix"],
    #     genecore_date_folder=config["genecore_date_folder"],
    # ),
    output_path="{data_location}/config/config_df.tsv".format(
        data_location=config["data_location"]
    ),
    check_sm_tag=False,
    bam=bam,
    genecore=config["genecore"],
)
# df_config_files = c.df_config_files
if config["genecore"] is True:
    d_master = c.d_master

# Read config file previously produced
df_config_files = c.df_config_files
df_config_files["Selected"] = True

# List of available samples
samples = list(sorted(list(df_config_files.Sample.unique().tolist())))


# scNOVA dedicated - to handle only selected cells labeled by the user /& ashleys
if config["scNOVA"] is True:
    l = list()
    # print(samples)
    for sample in samples:
        # Path of the labels file
        labels_path = "{folder}/{sample}/cell_selection/labels.tsv".format(
            folder=config["data_location"], sample=sample
        )

        assert os.path.isfile(
            labels_path
        ), "Ashleys labels were not computed yet, use first ashleys mode to perform cell selection"

        # print(labels_path)
        if os.path.exists(labels_path):
            # Read df
            tmp_df_labels_selected = pd.read_csv(labels_path, sep="\t")[
                ["cell", "prediction"]
            ]
            # Reformat to match #df_config_files
            tmp_df_labels_selected = tmp_df_labels_selected.rename(
                {"cell": "Cell", "prediction": "Selected"}, axis=1
            )
            tmp_df_labels_selected["Cell"] = tmp_df_labels_selected["Cell"].str.replace(
                ".sort.mdup.bam", "", regex=False
            )
            tmp_df_labels_selected["Selected"] = tmp_df_labels_selected[
                "Selected"
            ].astype(bool)
            # print(tmp_df_labels_selected)
            # print(df_config_files)
            # Merge dfs
            tmp_merge_df = pd.merge(
                tmp_df_labels_selected,
                df_config_files.drop(["Selected"], axis=1),
                on=["Cell"],
            )
            # Handle use-case if df don't have the same shapes
            if (
                tmp_merge_df.shape[0]
                < df_config_files.loc[df_config_files["Sample"] == sample].shape[0]
            ):
                print("WARNING: shape error when merging labels TSV & config TSV")
                tmp_merge_df = df_config_files.loc[
                    df_config_files["Sample"] == sample, ["Cell"]
                ]
                tmp_merge_df["Selected"] = True
            l.append(tmp_merge_df)
    # print(l)
    # Concat df to create a new one
    df_config_files_with_labels = pd.concat(l).reset_index(drop=True)
    df_config_files_with_labels.to_csv(
        "{data_location}/config/config_df_scNOVA.tsv".format(
            data_location=config["data_location"]
        )
    )

    bam_per_sample_selected = (
        df_config_files_with_labels.loc[df_config_files_with_labels["Selected"] == True]
        .groupby("Sample")["Cell"]
        .unique()
        .apply(list)
        .to_dict()
    )


# Creation of dicts to be used in the rules
dict_cells_nb_per_sample = (
    df_config_files.loc[df_config_files["Selected"] == True]
    .groupby("Sample")["Cell"]
    .nunique()
    .to_dict()
)

allbams_per_sample = df_config_files.groupby("Sample")["Cell"].apply(list).to_dict()
cell_per_sample = (
    df_config_files.loc[df_config_files["Selected"] == True]
    .groupby("Sample")["Cell"]
    .unique()
    .apply(list)
    .to_dict()
)
bam_per_sample_local = (
    df_config_files.loc[df_config_files["Selected"] == True]
    .groupby("Sample")["Cell"]
    .unique()
    .apply(list)
    .to_dict()
)

bam_per_sample = (
    df_config_files.loc[df_config_files["Selected"] == True]
    .groupby("Sample")["Cell"]
    .unique()
    .apply(list)
    .to_dict()
)


plottype_counts = (
    config["plottype_counts"]
    if config["multistep_normalisation"] is True
    else config["plottype_counts"][0]
)

if config["scNOVA"] is True:
    clones = collections.defaultdict(list)
    for sample in samples:
        subclonality_file = pd.read_csv(
            "{}/{}/scNOVA_input_user/input_subclonality.txt".format(
                config["data_location"], sample
            ),
            sep="\t",
        )
        clones[sample] = list(sorted(subclonality_file.Subclonality.unique().tolist()))
# print(clones)


def get_mem_mb(wildcards, attempt):
    mem_avail = [2, 4, 8, 16, 64]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    mem_avail = [8, 16, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000


def onsuccess_fct(log):
    config_metadata = config_definitions = yaml.safe_load(
        open(configfile_location.replace("config.yaml", "config_metadata.yaml"), "r")
    )
    log_path_new = make_log_useful.make_log_useful(
        log, "SUCCESS", config, config_metadata
    )
    shell(
        'mail -s "[smk-wf-catalog/mosaicatcher-pipeline] v{} - [{}--{}] - SUCCESS" {} < {}'.format(
            config["version"],
            config["data_location"].split("/")[-1],
            ";".join(samples),
            config["email"],
            log_path_new,
        )
    )


def onerror_fct(log):
    config_metadata = config_definitions = yaml.safe_load(
        open(configfile_location.replace("config.yaml", "config_metadata.yaml"), "r")
    )
    log_path_new = make_log_useful.make_log_useful(
        log, "ERROR", config, config_metadata
    )
    shell(
        'mail -s "[smk-wf-catalog/mosaicatcher-pipeline] v{} - [{}--{}] - ERROR" {} < {}'.format(
            config["version"],
            config["data_location"].split("/")[-1],
            ";".join(samples),
            config["email"],
            log_path_new,
        )
    )


def get_scnova_final_output(wildcards):
    # WARNING

    # subclonality_file = pd.read_csv(config["scnova_subclonality"], sep="\t")
    # clones = ["clone1", "clone2"]
    # clones = config["subclonality_list"]

    # SAMPLE_NAME = "TALL03-DEA5"
    # cell_per_sample[wildcards.sample], = glob_wildcards("input_bam/{cell}.bam")
    # abbreviate_names = False

    l = [
        expand(
            "{folder}/{sample}/scNOVA_result_plots/Result_scNOVA_plots_{sample}_alternative_PLSDA.pdf",
            folder=config["data_location"],
            sample=wildcards.sample,
        ),
        expand(
            "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort_num_sort_for_chromVAR.txt",
            folder=config["data_location"],
            sample=wildcards.sample,
        ),
        expand(
            "{folder}/{sample}/scNOVA_result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt",
            folder=config["data_location"],
            sample=wildcards.sample,
        ),
        expand(
            "{folder}/{sample}/scNOVA_result_haplo/Deeptool_Genebody_H1H2_sort.txt",
            folder=config["data_location"],
            sample=wildcards.sample,
        ),
    ]
    l = [sub_e for e in l for sub_e in e]
    return l


def get_final_output():
    """
    Input function of the pipeline, will retrieve all 'end' outputs
    """
    final_list = list()

    final_list.extend(
        expand(
            "{folder}/{sample}/plots/final_results/{sample}.txt",
            folder=config["data_location"],
            sample=samples,
        )
    )
    # print(config["scNOVA"])
    if config["scNOVA"] is True:
        # print("TOTO")
        final_list.extend(get_final_output_scnova())
    # from pprint import pprint
    # pprint(final_list)

    return final_list


def get_final_output_scnova():
    """
    Input function of the pipeline, will retrieve all 'end' outputs
    """
    final_list = list()

    final_list.extend(
        expand(
            "{folder}/{sample}/plots/final_results/scNOVA_{sample}.txt",
            folder=config["data_location"],
            sample=samples,
        )
    )

    return final_list


def get_all_plots(wildcards):
    """
    Function to retrieve all the plots/stats/outputs produced during the pipeline
    """

    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )

    dict_cells_nb_per_sample = {k: len(v) for k, v in cell_per_sample.items()}
    samples = list(dict_cells_nb_per_sample.keys())

    # QC Counts section
    # Create a tmp dictionnary and a corresponding PDF page for each of the cell of the run

    tmp_dict = {
        s: {i + 1: c for i, c in enumerate(cell_list)}
        for s, cell_list in cell_per_sample.items()
    }
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"

    l_outputs = list()

    tmp_l_divide_count_plots = [
        expand(
            "{folder}/{sample}/plots/counts/CountComplete.{plottype}.pdf",
            folder=config["data_location"],
            sample=sample,
            # plottype=["raw"],
            plottype=plottype_counts,
        )
        for sample in samples
    ]
    l_outputs.extend([sub_e for e in tmp_l_divide_count_plots for sub_e in e])

    if config["split_qc_plot"] is True:
        # print("OK")

        tmp_l_divide = [
            expand(
                "{folder}/{sample}/plots/counts_{plottype}/{cell}.{i}.pdf",
                folder=config["data_location"],
                sample=sample,
                plottype=["raw"],
                # plottype=plottype_counts,
                cell=tmp_dict[sample][i],
                i=i,
            )
            for sample in samples
            for i in range(dict_cells_nb_per_sample[sample] + 1)
        ]
        l_outputs.extend([sub_e for e in tmp_l_divide for sub_e in e])

    if config["arbigent"] is True:
        l_outputs.extend(
            expand(
                "{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/qc/lineplot_gts.pdf",
                folder=config["data_location"],
                sample=samples,
            )
        )

    else:
        # SV_consistency section

        l_outputs.extend(
            [
                sub_e
                for e in [
                    expand(
                        "{folder}/{sample}/plots/sv_consistency/{method}_filter{filter}.consistency-barplot-{plottype}.pdf",
                        folder=config["data_location"],
                        sample=wildcards.sample,
                        method=method,
                        plottype=config["plottype_consistency"],
                        filter=config["methods"][method]["filter"],
                    )
                    for method in config["methods"]
                ]
                for sub_e in e
            ]
        )

        l_outputs.extend(
            [
                sub_e
                for e in [
                    expand(
                        "{folder}/{sample}/plots/sv_clustering_dev/{method}-filter{filter}-{plottype}.pdf",
                        folder=config["data_location"],
                        sample=wildcards.sample,
                        method=method,
                        plottype=config["plottype_clustering"],
                        # plottype=config["plottype_clustering"],
                        filter=config["methods"][method]["filter"],
                    )
                    for method in config["methods"]
                ]
                for sub_e in e
            ]
        )

        l_outputs.extend(
            [
                sub_e
                for e in [
                    expand(
                        "{folder}/{sample}/plots/sv_calls_dev/{method}_filter{filter}/{chrom}.pdf",
                        folder=config["data_location"],
                        sample=wildcards.sample,
                        method=method,
                        chrom=config["chromosomes"],
                        filter=config["methods"][method]["filter"],
                    )
                    for method in config["methods"]
                ]
                for sub_e in e
            ]
        )

        # Complex section

        l_outputs.extend(
            [
                sub_e
                for e in [
                    expand(
                        "{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv",
                        folder=config["data_location"],
                        sample=wildcards.sample,
                        method=method,
                        filter=config["methods"][method]["filter"],
                    )
                    for method in config["methods"]
                ]
                for sub_e in e
            ]
        ),

        # Ploidy section
        l_outputs.extend(
            expand(
                "{folder}/{sample}/plots/ploidy/{sample}.pdf",
                folder=config["data_location"],
                sample=wildcards.sample,
            ),
        )

        # scTRIP multiplot

        if config["scTRIP_multiplot"] == True:
            l_outputs.extend(
                expand(
                    "{folder}/{sample}/plots/scTRIP_multiplot_aggr.ok",
                    folder=config["data_location"],
                    sample=wildcards.sample,
                )
            )

        # UCSC + IGV

        if config["genome_browsing_files_generation"] == True:
            l_outputs.extend(
                expand(
                    "{folder}/{sample}/plots/IGV/{sample}_IGV_session.xml",
                    folder=config["data_location"],
                    sample=wildcards.sample,
                ),
            )
            l_outputs.extend(
                expand(
                    "{folder}/{sample}/plots/UCSC/{sample}.bedUCSC.gz",
                    folder=config["data_location"],
                    sample=wildcards.sample,
                ),
            )
            l_outputs.extend(
                expand(
                    "{folder}/{sample}/plots/JBROWSE/{sample}.ok",
                    folder=config["data_location"],
                    sample=wildcards.sample,
                ),
            )

        # Stats section

        l_outputs.extend(
            expand(
                "{folder}/{sample}/stats/stats-merged.html",
                folder=config["data_location"],
                sample=wildcards.sample,
            ),
        )

        # Config section

        l_outputs.extend(
            expand(
                "{folder}/{sample}/config/config.yaml",
                folder=config["data_location"],
                sample=wildcards.sample,
            ),
        )

        # Run summary section

        # l_outputs.extend(
        #     expand(
        #         "{folder}/{sample}/config/run_summary.txt",
        #         folder=config["data_location"],
        #         sample=wildcards.sample,
        #     ),
        # )

    # from pprint import pprint
    # pprint(l_outputs)
    return l_outputs
