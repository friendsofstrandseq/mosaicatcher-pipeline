import pandas as pd
from scripts.utils import handle_input, make_log_useful, pipeline_aesthetic_start
import os

# Solve LC_CTYPE issue

os.environ["LC_CTYPE"] = "C"


envvars:
    "LC_CTYPE",


# Start with aesthetic pipeline config presentation
# onstart:
#     pipeline_aesthetic_start.pipeline_aesthetic_start(config)


# Configure if handle_input needs to be based on bam or fastq
bam = True if config["ashleys_pipeline"] is False else False


# Create configuration file with samples
c = handle_input.HandleInput(
    input_path=config["data_location"],
    output_path="{data_location}/config/config_df.tsv".format(
        data_location=config["data_location"]
    ),
    check_sm_tag=False,
    bam=bam,
)

# Read config file previously produced
df_config_files = c.df_config_files
df_config_files["Selected"] = True

# List of available samples
samples = list(sorted(list(df_config_files.Sample.unique().tolist())))


# List of assertions to verify
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
    mem_avail = [2, 4, 8, 16, 64, 128, 256]
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

    # l_outputs.extend(
    #     [
    #         sub_e
    #         for e in [
    #             expand(
    #                 "{folder}/{sample}/plots/sv_clustering/{method}-filter{filter}-{plottype}.pdf",
    #                 folder=config["data_location"],
    #                 sample=samples,
    #                 method=method,
    #                 plottype=config["plottype_clustering"],
    #                 filter=config["methods"][method]["filter"],
    #             )
    #             for method in config["methods"]
    #         ]
    #         for sub_e in e
    #     ]
    # )

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
