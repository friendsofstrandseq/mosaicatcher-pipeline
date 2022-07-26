# from workflow.scripts.utils.utils import get_mem_mb

import math

# import pandas as pd
# config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# cell_per_sample = config_df.loc[config_df["Selected"] == True].groupby("Sample")["Cell"].apply(list).to_dict()

################################################################################
# Joint Segmentation                                                                 #
################################################################################


# CHECKME : @Marco mention on Gitlab
# CHECKME : parameters
rule segmentation:
    """
    rule fct: Identify breakpoints of futur SV based on normalized read counts
    input: mosaic [normalized] counts
    output: Segmentation tab file 
    """
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        seg_filter="{output_folder}/log/segmentation/{sample}/segmentation-per-cell.ok",

    output:
        "{output_folder}/segmentation/{sample}/{sample}.txt.fixme",
    log:
        "{output_folder}/log/segmentation/{sample}/{sample}.log",
    params:
        min_num_segs=lambda wc: math.ceil(200000 / float(config["window"])),  # bins to represent 200 kb
    # container:
    #     "library://weber8thomas/remote-build/mosaic:0.3"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input.counts} > {log} 2>&1
        """


# FIXME: This is a workaround because latest versions of "mosaic segment" don't compute the "bps" column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    """
    rule fct:
    input:
    output:
    """
    input:
        ancient("{output_folder}/segmentation/{sample}/{sample}.txt.fixme"),
    output:
        "{output_folder}/segmentation/{sample}/{sample}.txt",
    log:
        "{output_folder}/log/segmentation/{sample}/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    params:
        script="workflow/scripts/segmentation_scripts/fix_segmentation.awk",
        window=config["window"],
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v name={wildcards.sample} -v window={params.window} -f {params.script} {input} > {output}
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
        "{output_folder}/counts/{sample}/counts-per-cell/{cell}.txt.gz",
    output:
        "{output_folder}/segmentation/{sample}/segmentation-per-cell/{cell}.txt",
    log:
        "{output_folder}/log/segmentation/{sample}/segmentation-per-cell/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    params:
        # mc_command = config["mosaicatcher"],
        min_num_segs=lambda wc: math.ceil(200000 / float(config["window"])),  # bins to represent 200 kb
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """

# rule aggregate_cells_merge:
#     input:
#         aggregate_cells_segmentation,
#     output:
#         "{output_folder}/log/segmentation/{sample}/segmentation-per-cell.ok"
#     shell:
#         "touch {output}"


rule segmentation_selection:
    """
    rule fct:
    input: mosaic read counts (txt.gz) & stats info (.info) + joint & sc segmentation 
    output: initial_strand_state used for the following by strandphaser
    """
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        jointseg="{output_folder}/segmentation/{sample}/{sample}.txt",
        singleseg=aggregate_cells_segmentation,
        info="{output_folder}/counts/{sample}/{sample}.info",
        seg_filter="{output_folder}/log/segmentation/{sample}/segmentation-per-cell.ok",
    output:
        jointseg="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
        singleseg="{output_folder}/segmentation/{sample}/Selection_singleseg.txt",
        strand_states="{output_folder}/segmentation/{sample}/Selection_initial_strand_state",
    log:
        "{output_folder}/log/segmentation/segmentation_selection/{sample}.log",
    params:
        cellnames=lambda wc: ",".join([cell for cell in cell_per_sample[wc.sample] if cell in [e.split(snakemake.config["abs_path"])[-1].split(".")[0] for e in aggregate_cells_segmentation(wc)]]),
        sce_min_distance=config["sce_min_distance"],
        additional_sce_cutoff=config["additional_sce_cutoff"],
        min_diff_jointseg=config["min_diff_jointseg"],
        min_diff_singleseg=config["min_diff_singleseg"],
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        python workflow/scripts/segmentation_scripts/detect_strand_states.py \
            --sce_min_distance {params.sce_min_distance} \
            --sce_add_cutoff {params.additional_sce_cutoff} \
            --min_diff_jointseg {params.min_diff_jointseg} \
            --min_diff_singleseg {params.min_diff_singleseg} \
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
