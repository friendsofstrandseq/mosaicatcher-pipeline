print("###################################")
print("# MOSAICATCHER SNAKEMAKE PIPELINE #")
print("###################################")

import math
from collections import defaultdict
import pandas as pd
import os, sys
from pprint import pprint
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("7.4.1")


# TODO : check config/config CLI ... 
# configfile: "../config/config.yaml"
configfile: "config/config.yaml"
# print(config)

mode_selected = config["mode"].lower()
correct_modes = ["count", "segmentation", "mosaiclassifier", "download_data"]
assert mode_selected in correct_modes, "Wrong mode selected : {}\nFollowing list of modes are available : {}".format(config["mode"], ", ".join(correct_modes))

plot_option_selected = config["plot"]
assert type(plot_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 

dl_bam_example_option_selected = config["dl_bam_example"]
assert type(dl_bam_example_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 

dl_external_files_option_selected = config["dl_external_files"]
assert type(dl_external_files_option_selected) is bool, "Wrong plot option selected : {}\nPlease enter a valid value (True / False)".format(config["plot"]) 


print("Mode selected : {}".format(mode_selected))
print("Plots output enabled : {}".format(config["plot"]))


print('Input folder selected : {}'.format(config['input_bam_location']))
print('Output folder selected : {}'.format(config['output_location']))

if config["output_location"].endswith("/") is False:
    config["output_location"] +=  "/"
    # print("/ was missing at end of output folder, continuing ...")

if config["input_bam_location"].endswith("/") is False:
    config["input_bam_location"] +=  "/"
    # print("/ was missing at end of input folder, continuing ...")



report: "workflow/report/workflow.rst"


# TODO: commit only in another branch to test
# # this container defines the underlying OS for each job when using the workflow
# # with --use-conda --use-singularity
# container: "docker://condaforge/mambaforge:4.12.0-0"

if config["containerized"] is True:
    containerized: "docker://weber8thomas/mosaicatcher-pipeline:1.5.1"

# include: "rules/input_check.smk"
# rule all:
#     input: 
#         rules.check_bam_input.output


if config["mode"] != "download_data":  
    if os.path.isfile(config["output_location"] + "config/config_df.tsv") is False:
        from workflow.scripts.utils import handle_input
        c = handle_input.HandleInput(
            input_path=config["input_bam_location"], 
            output_path=config["output_location"] + "config/config_df.tsv", 
            check_sm_tag=config["check_sm_tag"]
            )
        df_config_files = c.df_config_files
    else:
        df_config_files = pd.read_csv(config["output_location"] + "config/config_df.tsv", sep="\t")
    # print(df_config_files)
    # exit()




    ######################################################################

    # TODO : move to somewhere else

    methods = [
        "simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE",
        "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE",
    ]

    ######################################################################

    # IMPORT SMK RULES

    include: "workflow/rules/input_check.smk"
    include: "workflow/rules/setup.smk"
    include: "workflow/rules/count.smk"
    include: "workflow/rules/segmentation.smk"
    include: "workflow/rules/plots.smk"
    include: "workflow/rules/regenotyping.smk"
    include: "workflow/rules/utils.smk"
    include: "workflow/rules/strandphaser.smk"
    include: "workflow/rules/haplotagging.smk"
    include: "workflow/rules/mosaiclassifier.smk"
    include: "workflow/rules/postprocessing.smk"
    include: "workflow/rules/stats.smk"


include: "workflow/rules/examples.smk"

# Global wildcard constraints for consistent naming
# wildcard_constraints:



# IF PLOT OPTION ENABLED, BUILD TMP DICT TO CALL OUTPUT 
if plot_option_selected == True:
    dict_cells_nb_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["Cell"].nunique().to_dict()
    tmp_dict = df_config_files.loc[df_config_files["Selected"] == True, ["Sample", "Cell"]].groupby("Sample")["Cell"].apply(lambda r: sorted(list(r))).to_dict()
    tmp_dict = {s:{i+1:c for i,c in enumerate(cell_list)} for s,cell_list in tmp_dict.items()}
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"









######################
# MODES OF EXECUTION #
######################


# MODE MOSAIC COUNT
if mode_selected == "count":
    if plot_option_selected == True:
        rule all:
            input:
                [config["output_location"] + "plots/{}/counts/{}.{}.pdf".format(sample, tmp_dict[sample][i], i) for sample in samples for i in range(dict_cells_nb_per_sample[sample] + 1)],

    elif plot_option_selected == False:
        rule all:
            input:
                [config["output_location"] +  "counts/{}/{}.txt.fixme.gz".format(sample, sample) for sample in samples]

# MODE MOSAIC SEGMENTATION
elif mode_selected == "segmentation":

    if plot_option_selected == True:
        rule all:
            input:
                [config["output_location"] + "plots/{}/counts/{}.{}.pdf".format(sample, tmp_dict[sample][i], i) for sample in samples for i in range(dict_cells_nb_per_sample[sample] + 1)],
                [config["output_location"] + "segmentation/{}/Selection_initial_strand_state".format(sample) for sample in samples]
    
    elif plot_option_selected == False:
        rule all:
            input:
                [config["output_location"] + "segmentation/{}/Selection_initial_strand_state".format(sample) for sample in samples]

# MODE MOSAIC CLASSIFIER
elif mode_selected == "mosaiclassifier":
    if plot_option_selected == True:
        rule all:
            input:
                rules.cp_yaml_config.output,
                rules.install_rlib_strandphaser.output,
                [config["output_location"] + "plots/{}/counts/{}.{}.pdf".format(sample, tmp_dict[sample][i], i) for sample in samples for i in range(dict_cells_nb_per_sample[sample] + 1)],
                [config["output_location"] + "mosaiclassifier/sv_calls/{}/{}.tsv".format(sample, m) for sample in samples for m in methods],
                expand(config["output_location"] + "plots/{sample}/sv_calls/{method}.{chrom}.pdf", sample = samples, chrom = config["chromosomes"], method = methods),
                expand(config["output_location"] + "plots/{sample}/sv_consistency/{method}.consistency-barplot-{plottype}.pdf", sample = samples, method = methods, plottype = ["byaf","bypos"]),
                expand(config["output_location"] + "plots/{sample}/sv_clustering/{method}-{plottype}.pdf", sample = samples, method = methods, plottype = ["position", "chromosome"]),
                expand(config["output_location"] + "stats/{sample}/stats-merged.tsv", sample = samples),

    
    elif plot_option_selected == False:
        rule all:
            input:
                rules.install_rlib_strandphaser.output,
                [config["output_location"] + "mosaiclassifier/sv_calls/{}/{}.tsv".format(sample, m) for sample in samples for m in methods],
                [config["output_location"] + "mosaiclassifier/sv_calls/{}/{}.complex.tsv".format(sample, m) for sample in samples for m in methods],
                 expand(config["output_location"] + "stats/{sample}/stats-merged.tsv", sample = samples),



# # TEST MODE
elif mode_selected == "download_data":
    if dl_bam_example_option_selected is True and dl_external_files_option_selected is True:
        rule all:
            input:
                rules.dl_example_data.output,
                rules.dl_external_data.output,
                rules.dl_external_data_index.output
    
    if dl_bam_example_option_selected is True and dl_external_files_option_selected is False:
        rule all:
            input:
                rules.dl_example_data.output,
    
    if dl_bam_example_option_selected is False and dl_external_files_option_selected is True:
        rule all:
            input:
                rules.dl_external_data.output,
                rules.dl_external_data_index.output


if config["mail"]:
        
    onsuccess:
        print("Workflow finished, no error")
        shell('mail -s "Workflow finished, no error"' + config["mail"] + '< {log}')

    onerror:
        print("An error occurred")
        shell('mail -s "Workflow failed, an error occurred"' + config["mail"] + '< {log}')
