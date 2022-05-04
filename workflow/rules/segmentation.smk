import math
import pandas as pd
config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# print(config_df)
cell_per_sample = config_df.loc[config_df["all/selected"] == "selected"].groupby("Sample")["Cell"].apply(list).to_dict()

################################################################################
# Joint Segmentation                                                                 #
################################################################################


# CHECKME : @Marco mention on Gitlab
# CHECKME : parameters
# DOCME : check segmentation results to better understand
rule segmentation:
    """
    rule fct: Identify breakpoints of futur SV based on normalized read counts
    input: mosaic [normalized] counts
    output: Segmentation tab file 
    """
    input:
        config["output_location"] + "counts/{sample}/{sample}.txt.gz"
    output:
        config["output_location"] + "segmentation/{sample}/{sample}.txt.fixme"
    log:
        config["output_location"] + "log/segmentation/{sample}/{sample}.log"
    params:
        mc_command = config["mosaicatcher"],
        min_num_segs = lambda wc: math.ceil(200000 / float(config["window"]))  # bins to represent 200 kb
    container:
        "library://weber8thomas/remote-build/mosaic:0.3"
    shell:
        """
        /mosaicatcher/build/mosaic segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """
    # conda:
    #     "../envs/mc_base.yaml"
        # """
        # {params.mc_command} segment \
        # --remove-none \
        # --forbid-small-segments {params.min_num_segs} \
        # -M 50000000 \
        # -o {output} \
        # {input} 
        # """
        # {input} > {log} 2>&1

# FIXME: no difference observed before/after awk command


# FIXME: This is a workaround because latest versions of "mosaic segment" don't compute the "bps" column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["output_location"] + "segmentation/{sample}/{sample}.txt.fixme"
    output:
        config["output_location"] + "segmentation/{sample}/{sample}.txt"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v name={wildcards.sample} -v window={config[window]} -f scripts/segmentation_scripts/fix_segmentation.awk {input} > {output}
        """


################################################################################
# Single-Cell Segmentation                                                                 #
################################################################################


rule segment_one_cell:
    """
    rule fct: Same as `rule segmentation` : mosaic segment function but for individual cell
    input: mosaic count splitted by cell produced by `rule extract_single_cell_counts`
    output: Segmentation file for an individual cell
    """
    input:
        config["output_location"] + "counts/{sample}/counts-per-cell/{cell}.txt.gz"
    output:
        config["output_location"] + "segmentation/{sample}/segmentation-per-cell/{cell}.txt"
    log:
        config["output_location"] + "log/segmentation/{sample}/segmentation-per-cell/{cell}.log"
    container:
        "library://weber8thomas/remote-build/mosaic:0.3"
    params:
        # mc_command = config["mosaicatcher"],
        min_num_segs = lambda wc: math.ceil(200000 / float(config["window"])) # bins to represent 200 kb
    shell:
        """
        /mosaicatcher/build/mosaic segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """



rule segmentation_selection:
    """
    rule fct:
    input: mosaic read counts (txt.gz) & stats info (.info) + joint & sc segmentation 
    output: initial_strand_state used for the following by strandphaser
    """
    input:
        counts=config["output_location"] + "counts/{sample}/{sample}.txt.gz",
        jointseg=config["output_location"] + "segmentation/{sample}/{sample}.txt",
        singleseg=lambda wc: [config["output_location"] + "segmentation/{}/segmentation-per-cell/{}.txt".format(wc.sample, cell) for cell in cell_per_sample[wc.sample]],
        info=config["output_location"] + "counts/{sample}/{sample}.info",
    output:
        jointseg=config["output_location"] + "segmentation/{sample}/Selection_jointseg.txt",
        singleseg=config["output_location"] + "segmentation/{sample}/Selection_singleseg.txt",
        strand_states=config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state",
    log:
        config["output_location"] + "log/segmentation/segmentation_selection/{sample}.log"
    params:
        cellnames = lambda wc: ",".join(cell for cell in cell_per_sample[wc.sample]),
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        python scripts/segmentation_scripts/detect_strand_states.py \
            --sce_min_distance {config[sce_min_distance]} \
            --sce_add_cutoff {config[additional_sce_cutoff]} \
            --min_diff_jointseg {config[min_diff_jointseg]} \
            --min_diff_singleseg {config[min_diff_singleseg]} \
            --output_jointseg {output.jointseg} \
            --output_singleseg {output.singleseg} \
            --output_strand_states {output.strand_states} \
            --samplename {wildcards.sample} \
            --cellnames {params.cellnames} \
            {input.info} \
            {input.counts} \
            {input.jointseg} \
            {input.singleseg} > {log} 2>&1
        """

