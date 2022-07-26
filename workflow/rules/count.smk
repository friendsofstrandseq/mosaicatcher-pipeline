# from workflow.scripts.utils.utils import get_mem_mb

# import pandas as pd
# config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# bam_per_sample_local = config_df.loc[config_df["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()

# bam_per_sample_local = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()

################################################################################
# Read counting                                                                #
################################################################################


# DOCME : mosaic count read orientation ?
rule mosaic_count:
    """
    rule fct: Call mosaic count C++ function to count reads in each BAM file according defined window
    input: For the moment, individual BAM file in the selected folder of the associated sample
    output: counts: read counts for the BAM file according defined window ; info file : summary statistics 
    """
    input:
        bam=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam",
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)]
            if wc.sample in bam_per_sample_local
            else "FOOBAR",
        ),
        bai=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)],
        )
        if wc.sample in bam_per_sample_local
        else "FOOBAR",
        excl="{output_folder}/config_output/{sample}/exclude_file",
    output:
        counts="{output_folder}/counts/{sample}/{sample}.txt.raw.gz",
        info="{output_folder}/counts/{sample}/{sample}.info_raw",
    log:
        "{output_folder}/log/counts/{sample}/mosaic_count.log",
    # container:
    #     "library://weber8thomas/remote-build/mosaic:0.3"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    params:
        window=config["window"],
    resources:
        mem_mb=get_mem_mb,
    shell:
        # /mosaicatcher/build/mosaic count \
        """
        mosaicatcher count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {params.window} \
            {input.bam} \
        > {log} 2>&1
        """

if config["ashleys_pipeline"] is False:
    rule blank_labels:
        output:
            touch("{input_folder}/{sample}/cell_selection/labels.tsv")
        log:
            "{input_folder}/log/{sample}/blank_labels/labels.log"

rule order_mosaic_count_output:
    input:
        raw_count="{output_folder}/counts/{sample}/{sample}.txt.raw.gz",
        labels=expand(
            "{input_folder}/{sample}/cell_selection/labels.tsv",
            input_folder=config["input_bam_location"],
            sample=samples,
        ),
    output:
        "{output_folder}/counts/{sample}/{sample}.txt.sort.gz",
    log:
        "{output_folder}/log/counts/{sample}/{sample}.log",
    run:
        df = pd.read_csv(input.raw_count, compression="gzip", sep="\t")
        df = df.sort_values(by=["sample", "cell", "chrom", "start"])
        df.to_csv(output[0], index=False, compression="gzip", sep="\t")



checkpoint filter_bad_cells_from_mosaic_count:
    input:
        info_raw="{output_folder}/counts/{sample}/{sample}.info_raw",
        counts_sort="{output_folder}/counts/{sample}/{sample}.txt.sort.gz",
        labels=expand(
            "{input_folder}/{sample}/cell_selection/labels.tsv",
            input_folder=config["input_bam_location"],
            sample=samples,
        ),    
    output:
        info="{output_folder}/counts/{sample}/{sample}.info",
        info_removed="{output_folder}/counts/{sample}/{sample}.info_rm",
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
    log:
        "{output_folder}/log/filter_bad_cells_from_mosaic_count/{sample}/{sample}.log",
    script:
        "../scripts/utils/filter_bad_cells.py"


# CHECKME : to keep or to improve ? @jeong @mc @kg
################################################################################
# Normalize counts                                                             #
################################################################################


# TODO : check if inversion file is corresponded to previously published
rule merge_blacklist_bins:
    """
    rule fct: Call Python script to merge HGVSC normalization defined file & inversion whitelist file
    input: norm: HGSVC predefined BED file by the group ; whitelist: whitelist inversion file predefined by the group
    """
    input:
        norm="utils/normalization/HGSVC.{window}.txt",
        whitelist="utils/normalization/inversion-whitelist.tsv",
    output:
        merged="{output_folder}/normalizations/HGSVC.{window}.merged.tsv",
    log:
        "{output_folder}/log/merge_blacklist_bins/{window}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        utils/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2>> {log}
        """


# FIXME : snakemake ambiguity with I/O paths
# CHECKME : Check R code for normalization
rule normalize_counts:
    """
    rule fct: Normalization of mosaic counts based on merged normalization file produced with a linear relation (count * scaling_factor)
    input: counts: counts file coming from `rule mosaic_count` ; norm: merged normalization file produced by `rule merge_blacklist_bins`
    output: normalized counts based predefined factors for each window
    """
    input:
        counts="{output_folder}/counts/{sample}/{window}.txt.gz",
        norm="{output_folder}/normalizations/HGSVC.{window}.merged.tsv",
    output:
        "{output_folder}/norm_counts/{sample}/{window}.txt.gz",
    log:
        "{output_folder}/log/normalize_counts/{sample}/{window}.log",
    conda:
        "../envs/rtools.yaml"
    shell:
        """
        Rscript utils/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
        """


# FIXME : cleaner way to symlink info files
rule link_normalized_info_file:
    """
    rule fct: Symlink info file ouput mosaic count to normalization count directory
    input: Global summary statistics produced by mosaic count
    output: symlink in norm_counts output directory
    """
    input:
        info="{output_folder}/counts/{sample}/{window}.info",
    output:
        info="{output_folder}/norm_counts/{sample}/{window}.info",
    log:
        "{output_folder}/log/norm_counts/{sample}/{window}.log",
    run:
        d = os.path.dirname(output.info)
        file = os.path.basename(output.info)
        shell("cd {d} && ln -s {input.info} {file}")


################################################################################
# Single-Cell Segmentation                                                                 #
################################################################################


rule extract_single_cell_counts:
    """
    rule fct: extract from count the rows coming from the given cell
    input: mosaic count output file for the sample according a given window 
    output: count per cell file for the sample according a given window
    """
    input:
        info="{output_folder}/counts/{sample}/{sample}.info",
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        # agg=aggregate_cells
    output:
        "{output_folder}/counts/{sample}/counts-per-cell/{cell}.txt.gz",
    log:
        "{output_folder}/log/counts/{sample}/counts-per-cell/{cell}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        zcat {input.counts} | awk -v name={wildcards.cell} '(NR==1) || $5 == name' | gzip > {output}
        """
