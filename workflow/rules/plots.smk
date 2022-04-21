import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import pandas as pd
config_df = pd.read_csv("config/config_df.tsv", sep="\t")
tmp_dict = config_df.loc[config_df["all/selected"] == "selected", ["Sample", "Cell"]].groupby("Sample")["Cell"].apply(lambda r: sorted(list(r))).to_dict()
tmp_dict = {s:{i+1:c for i,c in enumerate(cell_list)} for s,cell_list in tmp_dict.items()}
for s in tmp_dict.keys():
    tmp_dict[s][0] = "SummaryPage"

print(tmp_dict)
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
            counts = config["output_location"] + "counts/{sample}/{sample}.txt.gz",
            info   = config["output_location"] + "counts/{sample}/{sample}.info"
        output:
            config["output_location"] + "plots/{sample}/CountComplete.pdf"
            # report(
            #     config["output_location"] + "plots/{sample}/Count_complete.pdf",
            #     category="Mosaic counts raw",
            # )
        log:
            config["output_location"] + "log/plot_mosaic_counts/{sample}.log"
        conda:
            "../envs/rtools.yaml"
        shell:
            """
            Rscript scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """

    rule divide_pdf:
        input:
            config["output_location"] + "plots/{sample}/CountComplete.pdf"
        output:
            report(
                config["output_location"] + "plots/{sample}/{cell}_{i, \d+}.pdf",
                    caption="../report/mosaic_counts.rst",
                    category="Mosaic counts",
            )
        run:

            from PyPDF2 import PdfFileWriter, PdfFileReader
            inputpdf = PdfFileReader(input[0], "rb")
            cell_name = tmp_dict[wildcards.sample][int(wildcards.i)] 
            output = PdfFileWriter()
            output.addPage(inputpdf.getPage(int(wildcards.i)))
            tmp_output_path = os.path.dirname(input[0]) + "/{}_{}.pdf".format(cell_name, wildcards.i)
            with open(tmp_output_path, "wb") as outputStream:
                output.write(outputStream)




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
