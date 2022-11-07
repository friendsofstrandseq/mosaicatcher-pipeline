import pandas as pd

# from scripts.utils import handle_input, make_log_useful, pipeline_aesthetic_start
from scripts.utils import make_log_useful, pipeline_aesthetic_start
import os, sys

# Solve LC_CTYPE issue

os.environ["LC_CTYPE"] = "C"


envvars:
    "LC_CTYPE",


# Start with aesthetic pipeline config presentation
onstart:
    pipeline_aesthetic_start.pipeline_aesthetic_start(config)


# List of assertions to verify
if config["genecore"] is False:
    l_samples = os.listdir(config["data_location"])
    assert (
        "fastq" not in l_samples
    ), "fastq folder found in the {} data_location specified: please specify a parent folder".format(
        config["data_location"]
    )

dl_bam_example_option_selected = config["dl_bam_example"]
assert (
    type(dl_bam_example_option_selected) is bool
), "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(
    config["plot"]
)

dl_external_files_option_selected = config["dl_external_files"]
assert (
    type(dl_external_files_option_selected) is bool
), "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(
    config["plot"]
)

if config["ashleys_pipeline"] is True:
    assert (
        config["ashleys_pipeline"] != config["input_bam_legacy"]
    ), "ashleys_pipeline and input_bam_legacy parameters cannot both be set to True"


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
                    lambda r: f"{config['genecore_prefix']}/{config['genecore_date_folder']}/{d_master[sample]['prefix']}lane1/{r.replace('.', '_')}_sequence.txt.gz"
                )
                df["Genecore_file"] = df["File"].apply(
                    lambda r: f"{d_master[sample]['prefix']}lane1{r.replace('.', '_')}"
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
                for f in os.listdir(
                    "{thisdir}/{sample}/{folder}/".format(
                        thisdir=thisdir, sample=sample, folder=folder
                    )
                )
                if f.endswith(ext)
            ]

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

            complete_df_list.append(df)

        # Concat dataframes for each sample & output
        complete_df = pd.concat(complete_df_list)
        complete_df = complete_df.sort_values(by=["Cell", "File"]).reset_index(
            drop=True
        )
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

c = HandleInput(
    input_path=config["data_location"],
    genecore_path="{genecore_prefix}/{genecore_date_folder}".format(
        genecore_prefix=config["genecore_prefix"],
        genecore_date_folder=config["genecore_date_folder"],
    ),
    output_path="{data_location}/config/config_df_ashleys.tsv".format(
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
    if config["GC_analysis"] is True
    else config["plottype_counts"][0]
)


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

    return final_list


def get_mem_mb(wildcards, attempt):
    mem_avail = [2, 4, 8, 16, 64]
    return mem_avail[attempt - 1] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    mem_avail = [8, 16, 64, 128, 256]
    return mem_avail[attempt - 1] * 1000


def onsuccess_fct(log):
    make_log_useful.make_log_useful(log, "SUCCESS", config)
    shell(
        'mail -s "[Snakemake] smk-wf-catalog/mosacaitcher-pipeline v{} - Run on {} - SUCCESS" {} < {{log}}'.format(
            config["version"], config["data_location"], config["email"]
        )
    )


def onerror_fct(log):
    make_log_useful.make_log_useful(log, "ERROR", config)
    shell(
        'mail -s "[Snakemake] smk-wf-catalog/mosacaitcher-pipeline v{} - Run on {} - ERRROR" {} < {{log}}'.format(
            config["version"], config["data_location"], config["email"]
        )
    )


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

    tmp_l_divide = [
        expand(
            "{folder}/{sample}/plots/counts_{plottype}/{cell}.{i}.pdf",
            folder=config["data_location"],
            sample=sample,
            plottype=plottype_counts,
            cell=tmp_dict[sample][i],
            i=i,
        )
        for sample in samples
        for i in range(dict_cells_nb_per_sample[sample] + 1)
    ]

    l_outputs.extend([sub_e for e in tmp_l_divide for sub_e in e])

    # SV_calls section

    l_outputs.extend(
        [
            sub_e
            for e in [
                expand(
                    "{folder}/{sample}/plots/sv_calls/{method}_filter{filter}/{chrom}.pdf",
                    folder=config["data_location"],
                    sample=samples,
                    method=method,
                    chrom=config["chromosomes"],
                    filter=config["methods"][method]["filter"],
                )
                for method in config["methods"]
            ]
            for sub_e in e
        ]
    )

    # SV_consistency section

    l_outputs.extend(
        [
            sub_e
            for e in [
                expand(
                    "{folder}/{sample}/plots/sv_consistency/{method}_filter{filter}.consistency-barplot-{plottype}.pdf",
                    folder=config["data_location"],
                    sample=samples,
                    method=method,
                    plottype=config["plottype_consistency"],
                    filter=config["methods"][method]["filter"],
                )
                for method in config["methods"]
            ]
            for sub_e in e
        ]
    )

    # SV_clustering section

    l_outputs.extend(
        [
            sub_e
            for e in [
                expand(
                    "{folder}/{sample}/plots/sv_clustering/{method}-filter{filter}-{plottype}.pdf",
                    folder=config["data_location"],
                    sample=samples,
                    method=method,
                    plottype=config["plottype_clustering"],
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
                    sample=samples,
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
            sample=samples,
        ),
    )

    # Stats section

    l_outputs.extend(
        expand(
            "{folder}/{sample}/stats/stats-merged.tsv",
            folder=config["data_location"],
            sample=samples,
        ),
    )

    # Run summary section

    l_outputs.extend(
        expand(
            "{folder}/config/{sample}/run_summary.txt",
            folder=config["data_location"],
            sample=samples,
        ),
    )
    return l_outputs
