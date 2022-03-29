import math
from collections import defaultdict

configfile: "Snake.config_embl.yaml"

import os, sys

# print(os.listdir(os.getcwd()))
# print(os.listdir("bam"))

# TODO I/O : Function to define inputs ; simplify list/dict system
# TODO Use remote file system to download example files

SAMPLE,BAM = glob_wildcards(config["input_bam_location"] + "{sample}/selected/{bam}.bam")

SAMPLES = sorted(set(SAMPLE))

CELL_PER_SAMPLE= defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)


for sample,bam in zip(SAMPLE,BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace(".sort.mdup",""))


ALLBAMS_PER_SAMPLE = defaultdict(list)
for sample in SAMPLES:
    ALLBAMS_PER_SAMPLE[sample] = glob_wildcards(config["input_bam_location"] + "{}/all/{{bam}}.bam".format(sample)).bam


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

print(BPDENS)
# Todo: specify an exact version of the singularity file!


print(SAMPLES)
print(CELL_PER_SAMPLE)
print(CELL_PER_SAMPLE.values())
print([sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e])
# print(expand([SAMPLES, [sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e]]))
print(expand(["{sample}/{cell}"], zip, sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e]))
# exit()
print(expand([config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz"], zip, sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e], window=[100000], ))

# rule all:
    # input:
        # expand(config["output_location"] + "counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "plots/{sample}/{window}.pdf", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "norm_counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
        # expand(config["output_location"] + "norm_counts/{sample}/{window}.info", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "segmentation/{sample}/{window}.txt", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "snv_genotyping/{sample}/{chrom}.vcf", sample=SAMPLES, window=[100000], chrom=config["chromosomes"])
        # expand(config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz", sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e], window=[100000], )
# FIXME : To solve : cell wildcard (dict type) comparatively to others that are list type


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
# DOCME : mosaic count read orientation ?
rule mosaic_count:
    """
    rule fct: Call mosaic count C++ function to count reads in each BAM file according defined window
    input: For the moment, individual BAM file in the selected folder of the associated sample
    output: counts: read counts for the BAM file according defined window ; info file : summary statistics 
    """
    input:
        bam = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        # excl = "log/exclude_file"
    output:
        counts = config["output_location"] + "counts/{sample}/{window}.txt.gz",
        info   = config["output_location"] + "counts/{sample}/{window}.info"
    log:
        config["output_location"] + "log/{sample}/mosaic_count.{window}.log"
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
    rule fct: Plot function of read counts for each bam file
    input: mosaic count outputs (counts & info)
    output: Generate figure based on couting results
    """
    input:
        counts = config["output_location"] + "counts/{sample}/{window}.txt.gz",
        info   = config["output_location"] + "counts/{sample}/{window}.info"
    output:
        config["output_location"] + "plots/{sample}/{window}.pdf"
    log:
        config["output_location"] + "log/plot_mosaic_counts/{sample}/{window}.log"
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
    rule fct: Normalization of mosaic counts based on merged normalization file produced 
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
        config["output_location"] + "counts/{sample}/{window}.txt.gz"
    output:
        config["output_location"] + "segmentation/{sample}/{window,\d+}.txt.fixme"
    log:
        config["output_location"] + "log/segmentation/{sample}/{window}.log"
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


# FIXME: This is a workaround because latest versions of "mosaic segment" don't compute the "bps" column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["output_location"] + "segmentation/{sample}/{window}.txt.fixme"
    output:
        config["output_location"] + "segmentation/{sample}/{window,\d+}.txt"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v name={wildcards.sample} -v window={wildcards.window} -f utils/command2.awk {input} > {output}
        """

# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["output_location"] + "segmentation/{sample}/{window}.txt"
    output:
        config["output_location"] + "segmentation2/{sample}/{window}.{bpdens,(many|medium|few)}.txt"
    log:
        config["output_location"] + "log/prepare_segments/{sample}/{window}.{bpdens}.log"
    params:
        quantile = lambda wc: config["bp_density"][wc.bpdens]
    script:
        "utils/helper.prepare_segments.R"

################################################################################
# Single-Cell Segmentation                                                                 #
################################################################################

# TODO : replace awk command with something else
rule extract_single_cell_counts:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["output_location"] + "counts/{sample}/{window}.txt.gz"
    output:
        config["output_location"] + "counts-per-cell/{sample}/{cell}/{window,[0-9]+}.txt.gz"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        zcat {input} | awk -v name={wildcards.cell} -f utils/command1.awk | gzip > {output}
        """

rule segment_one_cell:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz"
    output:
        config["output_location"] + "segmentation-per-cell/{sample}/{cell}/{window,\d+}.txt"
    log:
        config["output_location"] + "log/segmentation-per-cell/{sample}/{cell}/{window}.log"
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

rule segmentation_selection:
    input:
        counts=config["output_location"] + "counts/{sample}/{window}_{file_name}.txt.gz",
        jointseg=config["output_location"] + "segmentation/{sample}/{window}_{file_name}.txt",
        singleseg=lambda wc: [config["output_location"] + "segmentation-per-cell/{}/{}/{}_{}.txt".format(wc.sample, cell, wc.window, wc.file_name) for cell in CELL_PER_SAMPLE[wc.sample]],
        info=config["output_location"] + "counts/{sample}/{window}_{file_name}.info",
    output:
        jointseg=config["output_location"] + "segmentation2/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        singleseg=config["output_location"] + "segmentation-singlecell/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        strand_states=config["output_location"] + "strand_states/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}/intitial_strand_state",
    log:
        config["output_location"] + "log/segmentation_selection/{sample}/{window}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.log"
    params:
        cellnames = lambda wc: ",".join(cell for cell in CELL_PER_SAMPLE[wc.sample]),
        sce_min_distance = 500000,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        ./utils/detect_strand_states.py \
            --sce_min_distance {params.sce_min_distance} \
            --sce_add_cutoff {wildcards.additional_sce_cutoff}000000 \
            --min_diff_jointseg {wildcards.min_diff_jointseg} \
            --min_diff_singleseg {wildcards.min_diff_singleseg} \
            --output_jointseg {output.jointseg} \
            --output_singleseg {output.singleseg} \
            --output_strand_states {output.strand_states} \
            --samplename {wildcards.sample} \
            --cellnames {params.cellnames} \
            {input.info} \
            {input.counts} \
            {input.jointseg} \
            {input.singles
        """

################################################################################
# Call SNVs                                                                    #
################################################################################

# TODO : to move in utils category
rule mergeBams:
    """
    rule fct:
    input:
    output:
    """
    input:
        lambda wc: expand(config["input_bam_location"] + wc.sample + "/all/{bam}.bam", bam = ALLBAMS_PER_SAMPLE[wc.sample]) if wc.sample in ALLBAMS_PER_SAMPLE else "FOOBAR",
    output:
        config["output_location"] + "snv_calls/{sample}/merged.bam"
    log:
        config["output_location"] + "log/mergeBams/{sample}.log"
    threads:
        4
    shell:
        # FIXME : Samtools 1.10 from Conda env not working ; 1.9 from Seneca working > change it into conda env yml file
        config["samtools"] + " merge -@ {threads} {output} {input} 2>&1 > {log}"
        # "samtools" + " merge -@ {threads} {output} {input} 2>&1 > {log}"

# TODO : to move in utils category
rule index_bam:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.bam.log"
    shell:
        config["samtools"] + " index {input} 2> {log}"
        # "samtools" + " index {input} 2> {log}"


rule regenotype_SNVs:
    """
    rule fct:
    input:
    output:
    """
    input:
        bam   = config["output_location"] + "snv_calls/{sample}/merged.bam",
        bai   = config["output_location"] + "snv_calls/{sample}/merged.bam.bai",
        sites = config["snv_sites_to_genotype"],
    output:
        vcf = config["output_location"] + "snv_genotyping/{sample}/{chrom,chr[0-9A-Z]+}.vcf"
    log:
        config["output_location"] + "log/snv_genotyping/{sample}/{chrom}.log"
    params:
        fa = config["reference"],
        # bcftools = config["bcftools"]
    shell:
    # CHECKME : Samtools / BCFtools / freebayes path definition through conda env
    # CHECKME : interest of using -r parameters for freebayes => split by chroms
        """
        (freebayes \
            -f {params.fa} \
            -r {wildcards.chrom} \
            -@ {input.sites} \
            --only-use-input-alleles {input.bam} \
            --genotype-qualities \
        | bcftools view \
            --exclude-uncalled \
            --genotype het \
            --types snps \
            --include "QUAL>=10" - \
        > {output.vcf}) 2> {log}
        """
