import pandas as pd
from scripts.utils import handle_input, make_log_useful
import os

os.environ["LC_CTYPE"] = "C"


envvars:
    "LC_CTYPE",



onstart:
    sep = """------------------------------------------------------"""

    smk = """
                    _                        _        
    ___ _ __   __ _| | _____ _ __ ___   __ _| | _____ 
    / __| '_ \ / _` | |/ / _ \ '_ ` _ \ / _` | |/ / _ \\
    \__ \ | | | (_| |   <  __/ | | | | | (_| |   <  __/
    |___/_| |_|\__,_|_|\_\___|_| |_| |_|\__,_|_|\_\___|
    """

    wf_name = """                                                   
     __  __              _  ___      _      _            
    |  \/  |___ ___ __ _(_)/ __|__ _| |_ __| |_  ___ _ _ 
    | |\/| / _ (_-</ _` | | (__/ _` |  _/ _| ' \/ -_) '_|
    |_|  |_\___/__/\__,_|_|\___\__,_|\__\__|_||_\___|_|  
    """           
    wf_info = "smk-wf-catalog/mosacaitcher-pipeline v1.8.0"
    from colorama import Fore, Back, Style
    print(
        sep + Fore.GREEN + smk 
        )
    print(Style.RESET_ALL)
    print(Fore.YELLOW + wf_name)
    print(Style.RESET_ALL)
    print(Fore.MAGENTA + wf_info)
    print(Style.RESET_ALL)
    print(sep)

    # Input / Output
    print("\033[1m{}\033[0m".format("Input/Output options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Input folder selected", ": " + str(config["input_bam_location"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Output folder selected", ": " + str(config["output_location"])),
    ]
    [print(e) for e in l]

    print(Style.RESET_ALL)
    # Main options
    print("\033[1m{}\033[0m".format("Main options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("GC analysis", ": " + str(config["GC_analysis"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Read Counts normalization", ": " + str(config["normalized_counts"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Binning window size", ": " + str(config["window"])),
    ]
    [print(e) for e in l]
    
    print(Style.RESET_ALL)
    # Behavior options
    print("\033[1m{}\033[0m".format("Behavior options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Ashleys-QC preprocessing pipeline", ": " + str(config["ashleys_pipeline"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("BAM folder old format (all/selected)", ": " + str(config["input_old_behavior"])),
    ]
    [print(e) for e in l]
    
    print(Style.RESET_ALL)
    # Genome & chrom
    print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("List of chromosomes processed", ": " + ",".join(config["chromosomes"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Reference genome selected", ": " + str(config["reference"])),
        # f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Reference FASTA file", ": " + str(config["references_data"][config["reference"]]["reference_file_location"])),
    ]
    [print(e) for e in l]


if config["input_bam_location"] == config["output_location"]:
    sys.exit("input_bam_location & output_location must be different")


bam = True if config["ashleys_pipeline"] is False else False


# Create configuration file with samples
c = handle_input.HandleInput(
    input_path=config["input_bam_location"],
    output_path="{input_bam_location}/config/config_df.tsv".format(
        input_bam_location=config["input_bam_location"]
    ),
    check_sm_tag=False,
    bam=bam,
)

# Read config file previously produced
# df_config_files = pd.read_csv("{input_bam_location}/config/config_df.tsv".format(input_bam_location=config["input_bam_location"]), sep="\t")
df_config_files = c.df_config_files
df_config_files["Selected"] = True
# List of available samples
samples = list(sorted(list(df_config_files.Sample.unique().tolist())))
# print(df_config_files)


# mode_selected = config["mode"].lower()
# correct_modes = ["count", "segmentation", "mosaiclassifier", "download_data"]
# assert (
#     mode_selected in correct_modes
# ), "Wrong mode selected : {}\nFollowing list of modes are available : {}".format(
#     config["mode"], ", ".join(correct_modes)
# )

# plot_option_selected = config["plot"]
# assert (
#     type(plot_option_selected) is bool
# ), "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(
#     config["plot"]
# )

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
        config["ashleys_pipeline"] != config["input_old_behavior"]
    ), "ashleys_pipeline and input_old_behavior parameters cannot both be set to True"


dict_cells_nb_per_sample = (
    df_config_files.loc[df_config_files["Selected"] == True]
    .groupby("Sample")["Cell"]
    .nunique()
    .to_dict()
)

# samples = list(sorted(list(dict_cells_nb_per_sample.keys())))

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

samples_expand = [[k] * len(allbams_per_sample[k]) for k in allbams_per_sample.keys()]
samples_expand = [sub_e for e in samples_expand for sub_e in e]
cell_expand = [sub_e for e in list(allbams_per_sample.values()) for sub_e in e]

input_expand = [
    [config["input_bam_location"]] * len(allbams_per_sample[k])
    for k in allbams_per_sample.keys()
]
input_expand = [sub_e for e in input_expand for sub_e in e]

output_expand = [
    [config["output_location"]] * len(allbams_per_sample[k])
    for k in allbams_per_sample.keys()
]
output_expand = [sub_e for e in output_expand for sub_e in e]

plottype_counts = (
    config["plottype_counts"]
    if config["GC_analysis"] is True
    else config["plottype_counts"][0]
)


def get_final_output():

    final_list = list()

    final_list.extend(
        expand(
            "{output_folder}/plots/{sample}/final_results/{sample}.txt",
            output_folder=config["output_location"],
            sample=samples,
        )
    )

    return final_list


def get_mem_mb(wildcards, attempt):
    mem_avail = [2, 4, 8, 16, 64, 128, 256]
    # print(mem_avail[attempt-1] * 1000, attempt, mem_avail)
    return mem_avail[attempt - 1] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    mem_avail = [8, 16, 64, 128, 256]
    # print(mem_avail[attempt-1] * 1000, attempt, mem_avail)
    return mem_avail[attempt - 1] * 1000


def onsuccess_fct(wildcards):
    print("Workflow finished, no error")
    make_log_useful.make_log_useful(log, "SUCCESS", config["ouput_location", samples])

    shell('mail -s "[Snakemake] DGA - SUCCESS" {} < {{log}}'.format(config["email"]))


def onerror_fct(wildcards):
    print("An error occurred")
    make_log_useful.make_log_useful(log, "ERROR")

    shell('mail -s "[Snakemake] DGA - ERRROR" {} < {{log}}'.format(config["email"]))


def get_all_plots(wildcards):

    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, output_folder=config["output_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )

    dict_cells_nb_per_sample = df.groupby("sample")["cell"].nunique().to_dict()
    samples = list(dict_cells_nb_per_sample.keys())

    cell_list = df.cell.tolist()
    tmp_dict = (
        df[["sample", "cell"]]
        .groupby("sample")["cell"]
        .apply(lambda r: sorted(list(r)))
        .to_dict()
    )
    tmp_dict = {
        s: {i + 1: c for i, c in enumerate(cell_list)}
        for s, cell_list in tmp_dict.items()
    }
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"

    list_indiv_plots = list()

    tmp_l_divide = [
        expand(
            "{output_folder}/plots/{sample}/counts_{plottype}/{cell}.{i}.pdf",
            output_folder=config["output_location"],
            sample=sample,
            plottype=plottype_counts,
            cell=tmp_dict[sample][i],
            i=i,
        )
        for sample in samples
        for i in range(dict_cells_nb_per_sample[sample] + 1)
    ]

    list_indiv_plots.extend([sub_e for e in tmp_l_divide for sub_e in e])

    list_indiv_plots.extend(
        [
            sub_e
            for e in [
                expand(
                    "{output_folder}/plots/{sample}/sv_calls/{method}_filter{filter}/{chrom}.pdf",
                    output_folder=config["output_location"],
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
    list_indiv_plots.extend(
        [
            sub_e
            for e in [
                expand(
                    "{output_folder}/plots/{sample}/sv_consistency/{method}_filter{filter}.consistency-barplot-{plottype}.pdf",
                    output_folder=config["output_location"],
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
    list_indiv_plots.extend(
        [
            sub_e
            for e in [
                expand(
                    "{output_folder}/plots/{sample}/sv_clustering/{method}-filter{filter}-{plottype}.pdf",
                    output_folder=config["output_location"],
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
    list_indiv_plots.extend(
        [
            sub_e
            for e in [
                expand(
                    "{output_folder}/mosaiclassifier/complex/{sample}/{method}_filter{filter}.tsv",
                    output_folder=config["output_location"],
                    sample=samples,
                    method=method,
                    filter=config["methods"][method]["filter"],
                )
                for method in config["methods"]
            ]
            for sub_e in e
        ]
    ),
    list_indiv_plots.extend(
        expand(
            "{output_folder}/plots/{sample}/ploidy/{sample}.pdf",
            output_folder=config["output_location"],
            sample=samples,
        ),
    )
    # if config["GC_analysis"] is True:
    #     list_indiv_plots.extend(
    #         expand(
    #             "{output_folder}/plots/{sample}/alfred/gc_devi.png",
    #             output_folder=config["output_location"],
    #             sample=samples,
    #         ),
    #     )
    #     list_indiv_plots.extend(
    #         expand(
    #             "{output_folder}/plots/{sample}/alfred/gc_dist.png",
    #             output_folder=config["output_location"],
    #             sample=samples,
    #         ),
    #     )
    # list_indiv_plots.extend(
    #     expand(
    #         "{output_folder}/plots/{sample}/counts/CountComplete.{plottype_counts}.pdf",
    #         output_folder=config["output_location"],
    #         sample=samples,
    #         plottype_counts=config["plottype_counts"],
    #     ),
    # )
    list_indiv_plots.extend(
        expand(
            "{output_folder}/stats/{sample}/stats-merged.tsv",
            output_folder=config["output_location"],
            sample=samples,
        ),
    )
    list_indiv_plots.extend(
        expand(
            "{output_folder}/config/{sample}/run_summary.txt",
            output_folder=config["output_location"],
            sample=samples,
        ),
    )
    return list_indiv_plots
