from workflow.scripts.utils.utils import get_mem_mb 

import pandas as pd
config_df = pd.read_csv(config["output_location"] + "config/config_df.tsv", sep="\t")

pd.options.display.max_colwidth = 40
bam_per_sample_local = config_df.loc[config_df["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()

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
        bam = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam", bam = bam_per_sample_local[str(wc.sample)] if wc.sample in bam_per_sample_local else "FOOBAR"),
        # bai = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam.bai", bam = bam_per_sample_local[wc.sample]) if wc.sample in bam_per_sample_local else "FOOBAR",
        excl = ancient(config["output_location"] + "config/exclude_file"),
    output:
        counts = config["output_location"] + "counts/{sample}/{sample}.txt.fixme.gz",
        info   = config["output_location"] + "counts/{sample}/{sample}.info"
    log:
        config["output_location"] + "log/counts/{sample}/mosaic_count.log"
    container:
        "library://weber8thomas/remote-build/mosaic:0.3"
    params:
        window = config["window"]
    resources:
        mem_mb = get_mem_mb,
    shell:
        """
        /mosaicatcher/build/mosaic count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {params.window} \
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


