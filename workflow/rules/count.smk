import pandas as pd
config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# print(config_df)
bam_per_sample = config_df.loc[config_df["all/selected"] == "selected"].groupby("Sample")["File"].apply(list).to_dict()

################################################################################
# Read counting                                                                #
################################################################################

rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        "config/config_df.tsv"
    output:
        "config/exclude_file.txt"
    params:
        chroms = config["chromosomes"]
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/generate_exclude_file.py"

# TODO : find a solution to compile mosaic count automatically
# rule compile_mosaic:
#     """
#     rule fct:
#     input:
#     output:
#     """
#     input: 
#         "../src/"
#     output:
#         "../build/mosaic"
#     conda:
#         "../envs/mc_base.yaml"
#     shell:
#         """
#         mkdir -p ../build
#         cd ../build
#         cmake ../src
#         make
#         cd ../worfklow
#         """

# TODO : Simplify expand command 
# DOCME : mosaic count read orientation ?
rule mosaic_count:
    """
    rule fct: Call mosaic count C++ function to count reads in each BAM file according defined window
    input: For the moment, individual BAM file in the selected folder of the associated sample
    output: counts: read counts for the BAM file according defined window ; info file : summary statistics 
    """
    input:
        bam = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        # bai = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam.bai", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        excl = "config/exclude_file.txt",
        # mosaic = "../build/mosaic"
    output:
        counts = temp(config["output_location"] + "counts/{sample}/{sample}.txt.fixme.gz"),
        info   = config["output_location"] + "counts/{sample}/{sample}.info"
    log:
        config["output_location"] + "log/counts/{sample}/mosaic_count.log"
    conda:
        "../envs/mc_base.yaml"
    params:
        mc_command = config["mosaicatcher"],
        window = config["window"]
    shell:
        """
        {params.mc_command} count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -w {params.window} \
            -x {input.excl} \
            {input.bam} \
        > {log} 2>&1
        """


rule order_mosaic_count_output:
    input:
        config["output_location"] + "counts/{sample}/{sample}.txt.fixme.gz"
    output:
        config["output_location"] + "counts/{sample}/{sample}.txt.gz"

    run:
        df = pd.read_csv(input[0], compression='gzip', sep='\t')
        df = df.sort_values(by=["sample", "cell", "chrom", "start"])
        df.to_csv(output[0], index=False, compression="gzip", sep="\t")


# CHECKME : to keep or to improve ? @jeong @mc @kg
################################################################################
# Normalize counts                                                             #
################################################################################

# TODO : Reference blacklist BED file to retrieve easily on Git/Zenodo/remote system
# TODO : check if inversion file is corresponded to previously published 
rule merge_blacklist_bins:
    """
    rule fct: Call Python script to merge HGVSC normalization defined file & inversion whitelist file
    input: norm: HGSVC predefined BED file by the group ; whitelist: whitelist inversion file predefined by the group
    """
    input:
        norm = "utils/normalization/HGSVC.{window}.txt",
        whitelist = "utils/normalization/inversion-whitelist.tsv",
    output:
        merged = config["output_location"] + "normalizations/HGSVC.{window}.merged.tsv"
    log:
        config["output_location"] + "log/merge_blacklist_bins/{window}.log"
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
        counts = config["output_location"] + "counts/{sample}/{window}.txt.gz",
        norm   = config["output_location"] + "normalizations/HGSVC.{window}.merged.tsv",
    output:
        config["output_location"] + "norm_counts/{sample}/{window}.txt.gz"
    log:
        config["output_location"] + "log/normalize_counts/{sample}/{window}.log"
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
        info = config["output_location"] + "counts/{sample}/{window}.info"
    output:
        info = config["output_location"] + "norm_counts/{sample}/{window}.info"
    run:
        d = os.path.dirname(output.info)
        file = os.path.basename(output.info)
        shell("cd {d} && ln -s {input.info} {file}")



################################################################################
# Single-Cell Segmentation                                                                 #
################################################################################

# TODO : replace awk external file command with something else
rule extract_single_cell_counts:
    """
    rule fct: extract from count the rows coming from the given cell
    input: mosaic count output file for the sample according a given window 
    output: count per cell file for the sample according a given window
    """
    input:
        config["output_location"] + "counts/{sample}/{sample}.txt.gz"
    output:
        config["output_location"] + "counts/{sample}/counts-per-cell/{cell}.txt.gz"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        zcat {input} | awk -v name={wildcards.cell} '(NR==1) || $5 == name' | gzip > {output}
        """