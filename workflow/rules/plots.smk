import pandas as pd
config_df = pd.read_csv("config/config_df.tsv", sep="\t")
tmp_dict = config_df.loc[config_df["all/selected"] == "selected", ["Sample", "Cell"]].groupby("Sample")["Cell"].apply(lambda r: sorted(list(r))).to_dict()
tmp_dict = {s:{i+1:c for i,c in enumerate(cell_list)} for s,cell_list in tmp_dict.items()}
for s in tmp_dict.keys():
    tmp_dict[s][0] = "SummaryPage"
    
################################################################################
# Plots                                                                        #
################################################################################


# rule tmp_small_plots:
#     input:
#         counts = config["output_location"] + "counts/{sample}.txt.gz",
#         info   = config["output_location"] + "counts/{sample}.info"
#     output:
#         counts = config["output_location"] + "counts/{sample}.lite.txt.gz",
#         info   = config["output_location"] + "counts/{sample}.lite.info"
#     run:
#         df_counts = pd.read_csv(input[0], compression="gzip", sep="\t")
#         sample_cells = df_counts["cell"].unique().tolist()[:2]
#         df_counts_lite = df_counts.loc[df_counts["cell"].isin(sample_cells)]
#         df_counts_lite.to_csv(output[0], sep="\t", compression="gzip", index=False)

#         info_header = "".join([line for line in open(input[1], "r").readlines() if line[0] == "#"])
#         with open(output[1], "w") as w:
#             w.write(info_header)

#         df_info = pd.read_csv(input[1], sep="\t", skiprows=13)
#         df_info_lite = df_info.loc[df_info["cell"].isin(sample_cells)]
#         df_info_lite.to_csv(output[1], sep="\t", mode="a", index=False)


        

# FIXME : Missing plots in final PDF ; R script + inputs to check
# CHECKME : check if possible to switch from PDF to svg (or both) to produce lighter files
if config["plot"] is True:

    rule plot_mosaic_counts:
        """
        rule fct: Plot function of read counts for each bam file
        input: mosaic count outputs (counts & info)
        output: Generate figure based on couting results
        """
        input:
            counts = config["output_location"] + "counts/{sample}.txt.gz",
            info   = config["output_location"] + "counts/{sample}.info"
        output:
            config["output_location"] + "plots/{sample}/Count_complete.pdf"
            # report(
            #     config["output_location"] + "plots/{sample}/Count_complete.pdf",
            #     category="Mosaic counts raw",
            # )
        log:
            config["output_location"] + "log/plot_mosaic_counts/{sample}.log"
        conda:
            "../envs/plots.yaml"
        shell:
            """
            Rscript scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """

    ruleorder: divide_pdf_into_png > rename_png_cells

    rule divide_pdf_into_png:
        input:
            config["output_location"] + "plots/{sample}/Count_complete.pdf"
        output:
            config["output_location"] + "plots/{sample}/{i, \d+}.tmp.png"
        conda:
            "../envs/imagemagick.yaml"
        shell:
            'convert {input}"[{wildcards.i}]" {output}'


    rule rename_png_cells:
        input:
            config["output_location"] + "plots/{sample}/{i, \d+}.tmp.png"
        output:
            report(
                config["output_location"] + "plots/{sample}/{cell}_{i, \d+}.png",
                    caption="../report/mosaic_counts.rst",
                    category="Mosaic counts",
            )
        run:
            # import shutil
            input_name = input[0]
            cell_name = tmp_dict[wildcards.sample][int(wildcards.i)] 
            os.rename(input[0], os.path.dirname(input_name) + "/{}_{}.png".format(cell_name, wildcards.i))
            # shutil.copy(input[0], os.path.dirname(input_name) + "/{}_{}_{}.png".format(wildcards.sample, wildcards.i, cell_name))




    ## CHECKME : check halo ? 
    # rule generate_halo_json:
    #     input:
    #         counts = config["output_location"] + "counts/{sample}/{windows}.txt.gz",
    #     output:
    #         json = config["output_location"] + "halo/{sample}/{windows}.json.gz",
    #     log:
    #         config["output_location"] + "log/generate_halo_json/{sample}/{windows}.{windows}.log"
    #     shell:
    #         """
    #         PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
    #         (./utils/counts_to_json.py {input.counts} | gzip > {output.json}) 
    #         """
