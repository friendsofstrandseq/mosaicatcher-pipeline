import pandas as pd


# samples = (
#     pd.read_csv(config["samples"], sep="\t")
# )


mode_selected = config["mode"].lower()
correct_modes = ["count", "segmentation", "mosaiclassifier", "download_data"]
assert mode_selected in correct_modes, "Wrong mode selected : {}\nFollowing list of modes are available : {}".format(config["mode"], ", ".join(correct_modes))

plot_option_selected = config["plot"]
assert type(plot_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 

dl_bam_example_option_selected = config["dl_bam_example"]
assert type(dl_bam_example_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 

dl_external_files_option_selected = config["dl_external_files"]
assert type(dl_external_files_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 




df_config_files = pd.read_csv(config["samples"], sep="\t")
dict_cells_nb_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["Cell"].nunique().to_dict()
samples = list(sorted(list(dict_cells_nb_per_sample.keys())))
allbams_per_sample = df_config_files.groupby("Sample")["File"].apply(list).to_dict()
cell_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["Cell"].apply(list).to_dict()
bam_per_sample_local = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()
bam_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()

def get_final_output():


    final_list = list()
    # final_list.extend([config["output_location"] + "plots/{}/counts/{}.{}.pdf".format(sample, tmp_dict[sample][i], i) for sample in samples for i in range(dict_cells_nb_per_sample[sample] + 1)]),
    final_list.extend(expand("{output}/plots/{sample}/final_results/{sample}.txt", output=config['output_location'], sample = samples))
    # from pprint import pprint
    # pprint(final_list)
    return final_list



def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [1, 2, 4, 8, 16, 64, 128, 256 ]
    # print(mem_avail[attempt-1] * 1000, attempt, mem_avail)
    return mem_avail[attempt-1] * 1000


def get_indiv_plots_count():
    dict_cells_nb_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["Cell"].nunique().to_dict()
    samples = list(dict_cells_nb_per_sample.keys())
    tmp_dict = df_config_files.loc[df_config_files["Selected"] == True, ["Sample", "Cell"]].groupby("Sample")["Cell"].apply(lambda r: sorted(list(r))).to_dict()
    tmp_dict = {s:{i+1:c for i,c in enumerate(cell_list)} for s,cell_list in tmp_dict.items()}
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"
    list_indiv_plots = ["{}/plots/{}/counts/{}.{}.pdf".format(config['output_location'], sample, tmp_dict[sample][i], i) for sample in samples for i in range(dict_cells_nb_per_sample[sample] + 1)]

    list_indiv_plots.extend(expand("{output}/plots/{sample}/sv_calls/{method}.{chrom}.pdf", output=config['output_location'], sample = samples, chrom = config["chromosomes"], method = config["methods"]),)
    list_indiv_plots.extend(expand("{output}/plots/{sample}/sv_consistency/{method}.consistency-barplot-{plottype}.pdf", output=config['output_location'], sample = samples, method = config["methods"], plottype = config["plottype_consistency"]),)
    list_indiv_plots.extend(expand("{output}/plots/{sample}/sv_clustering/{method}-{plottype}.pdf", output=config['output_location'], sample = samples, method = config["methods"], plottype = config["plottype_clustering"]),)
    list_indiv_plots.extend(expand("{output}/stats/{sample}/stats-merged.tsv", output=config['output_location'], sample = samples),)
    list_indiv_plots.extend(expand("{output}/mosaiclassifier/sv_calls/{sample}/{method}.tsv", output=config['output_location'], sample = samples, method = config["methods"])),
    # print(list_indiv_plots)
    return list_indiv_plots

# def get_final_plots():
#     final_list = list()
#     return final_list