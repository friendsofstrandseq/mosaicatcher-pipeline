import math
from collections import defaultdict

configfile: "Snake.config_embl.yaml"

import os, sys

# print(os.listdir(os.getcwd()))
# print(os.listdir("bam"))

# TODO I/O : Function to define inputs ; simplify list/dict system
# TODO Use remote file system to download example files

SAMPLE,BAM = glob_wildcards(config['input_bam_location'] + "{sample}/selected/{bam}.bam")

SAMPLES = sorted(set(SAMPLE))

CELL_PER_SAMPLE= defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)


for sample,bam in zip(SAMPLE,BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace('.sort.mdup',''))


ALLBAMS_PER_SAMPLE = defaultdict(list)
for sample in SAMPLES:
    ALLBAMS_PER_SAMPLE[sample] = glob_wildcards(config['input_bam_location'] + "{}/all/{{bam}}.bam".format(sample)).bam


print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print("  {}:\t{} cells\t {} selected cells".format(s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])))


# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher


METHODS = [
    "simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE",
]

BPDENS = [
    "selected_j{}_s{}_scedist{}".format(joint, single, scedist) for joint in [0.1] for single in [0.5] for scedist in [20]
]

# Todo: specify an exact version of the singularity file!



rule all:
    input:
        # expand(config['output_location'] + "counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000])
        # expand(config['output_location'] + "plots/{sample}/{window}.pdf", sample=SAMPLES, window=[100000])
        # expand(config['output_location'] + "norm_counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
        # expand(config['output_location'] + "norm_counts/{sample}/{window}.info", sample=SAMPLES, window=[100000])
        expand(config['output_location'] + "segmentation/{sample}/{window}.txt", sample=SAMPLES, window=[100000])



################################################################################
# Read counting                                                                #
################################################################################

# CHECKME : exclude file rule useful ?
# rule generate_exclude_file_1:
#     output:
#         temp("log/exclude_file.temp")
#     input:
#         bam = expand("bam/{sample}/selected/{bam}.bam", sample = SAMPLES[0], bam = BAM_PER_SAMPLE[SAMPLES[0]][0])
#     log:
#         "log/generate_exclude_file_1.log"
#     params:
#         samtools = config["samtools"]
#     shell:
#         """
#         {params.samtools} view -H {input.bam} | awk "/^@SQ/" > {output} 
#         """

# rule generate_exclude_file_2:
#     output:
#         "log/exclude_file"
#     input:
#         "log/exclude_file.temp"
#     params:
#         chroms = config["chromosomes"]
#     run:
#         with open(input[0]) as f:
#             with open(output[0],"w") as out:
#                 for line in f:
#                     contig = line.strip().split()[1]
#                     contig = contig[3:]
#                     # if contig not in params.chroms:
#                         # print(contig, file = out)
# CHECKME : same as above for input ???
# TODO : Simplify expand command 
rule mosaic_count:
    """
    Call mosaic count C++ function to count reads in each BAM file according defined window
    """
    input:
        bam = lambda wc: expand(config['input_bam_location'] + wc.sample +  "/selected/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        # excl = "log/exclude_file"
    output:
        counts = config['output_location'] + "counts/{sample}/{window}.txt.gz",
        info   = config['output_location'] + "counts/{sample}/{window}.info"
    log:
        config['output_location'] + "log/{sample}/mosaic_count.{window}.log"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -w {wildcards.window} \
            {input.bam} 
        > {log} 2>&1
        """

# FIXME : Missing plots in final PDF ; R script + inputs to check
rule plot_mosaic_counts:
    """
    Plot function of read counts for each bam file
    """
    input:
        counts = config['output_location'] + "counts/{sample}/{window}.txt.gz",
        info   = config['output_location'] + "counts/{sample}/{window}.info"
    output:
        config['output_location'] + "plots/{sample}/{window}.pdf"
    log:
        config['output_location'] + "log/plot_mosaic_counts/{sample}/{window}.log"
    params:
        plot_command = "Rscript " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output} > {log} 2>&1
        """


################################################################################
# Normalize counts                                                             #
################################################################################

# TODO : Reference blacklist BED file to retrieve easily on Git/Zenodo/remote system
rule merge_blacklist_bins:
    """
    Call Python script to merge HGVSC normalization defined file & inversion whitelist file
    """
    input:
        norm = "utils/normalization/HGSVC.{window}.txt",
        whitelist = "utils/normalization/inversion-whitelist.tsv",
    output:
        merged = config['output_location'] + "normalizations/HGSVC.{window}.merged.tsv"
    log:
        config['output_location'] + "log/merge_blacklist_bins/{window}.log"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        utils/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2>> {log}
        """

# FIXME : snakemake ambiguity with I/O paths
# CHECKME : Check R code for normalization ; @Marco mention on Gitlab
rule normalize_counts:
    """
    Normalization of mosaic counts based on normalization file produced above
    """
    input:
        counts = config['output_location'] + "counts/{sample}/{window}.txt.gz",
        norm   = config['output_location'] + "normalizations/HGSVC.{window}.merged.tsv",
    output:
        config['output_location'] + "norm_counts/{sample}/{window}.txt.gz"
    log:
        config['output_location'] + "log/normalize_counts/{sample}/{window}.log"
    shell:
        """
        Rscript utils/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
        """

# FIXME : cleaner way to symlink info files
rule link_normalized_info_file:
    """
    Symlink info file ouput mosaic count to normalization count directory
    """
    input:
        info = config['output_location'] + "counts/{sample}/{window}.info"
    output:
        info = config['output_location'] + "norm_counts/{sample}/{window}.info"
    run:
        d = os.path.dirname(output.info)
        file = os.path.basename(output.info)
        shell("cd {d} && ln -s {input.info} {file}")



################################################################################
# Segmentation                                                                 #
################################################################################

# CHECKME : @Marco mention on Gitlab
rule segmentation:
    """
    rule fct:
    input:
    output:
    """
    input:
        config['output_location'] + "counts/{sample}/{window}.txt.gz"
    output:
        config['output_location'] + "segmentation/{sample}/{window,\d+}.txt"
    log:
        config['output_location'] + "log/segmentation/{sample}/{window}.log"
    params:
        mc_command = config["mosaicatcher"],
        min_num_segs = lambda wc: math.ceil(200000 / float(wc.window)) # bins to represent 200 kb
    shell:
        """
        {params.mc_command} segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """