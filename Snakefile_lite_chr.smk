import math
from collections import defaultdict

configfile: "Snake.config_embl_chr.yaml"
import pandas as pd
import os, sys
from pprint import pprint

# print(os.listdir(os.getcwd()))
# print(os.listdir("bam"))

# TODO I/O : Function to define inputs ; simplify list/dict system
# TODO Use remote file system to download example files

SAMPLE,BAM = glob_wildcards(config["input_bam_location_raw"] + "{sample}/selected/{bam}.bam")
# pprint(SAMPLE)
# pprint(BAM)
SAMPLES = sorted(set(SAMPLE))

CELL_PER_SAMPLE= defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)


for sample,bam in zip(SAMPLE,BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace(".sort.mdup",""))


ALLBAMS_PER_SAMPLE = defaultdict(list)
for sample in SAMPLES:
    ALLBAMS_PER_SAMPLE[sample] = glob_wildcards(config["input_bam_location_raw"] + "{}/all/{{bam}}.bam".format(sample)).bam


print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print("  {}:\t{} cells\t {} selected cells".format(s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])))

# pprint(BAM_PER_SAMPLE)
# pprint(CELL_PER_SAMPLE)

# FIXME : tmp solution to remove bad cells => need to fix this with combination of ASHLEYS ?

exclude_list = ['BM510x3PE20490']

BAM_PER_SAMPLE = {k:sorted([e  for e in v if e.split('.')[0] not in exclude_list]) for k,v in BAM_PER_SAMPLE.items()}
CELL_PER_SAMPLE = {k:sorted([e for e in v if e not in exclude_list]) for k,v in CELL_PER_SAMPLE.items()}


# pprint(BAM_PER_SAMPLE)
# pprint(CELL_PER_SAMPLE)

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

# FIXME : move to yaml/json settings or to something else
BPDENS = [
    "selected_j{}_s{}_scedist{}".format(joint, single, scedist) for joint in [0.1] for single in [0.5] for scedist in [20]
]

# print(BPDENS)
# # Todo: specify an exact version of the singularity file!


# print(SAMPLES)
# print(CELL_PER_SAMPLE)
# print(CELL_PER_SAMPLE.values())
# print([sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e])
# # print(expand([SAMPLES, [sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e]]))
# print(expand(["{sample}/{cell}"], zip, sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e]))
# # exit()
# print(expand([config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz"], zip, sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e], window=[100000], ))


final_list = [config['input_bam_location_raw'] + "{}/{}.bam".format(key, nested_key) for key in BAM_PER_SAMPLE for nested_key in BAM_PER_SAMPLE[key] ]



rule all:
    input:
        expand(config["output_location"] + "counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
        # expand(config["output_location"] + "plots/{sample}/{window}.pdf", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "norm_counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
        # expand(config["output_location"] + "norm_counts/{sample}/{window}.info", sample=SAMPLES, window=[100000])
        # expand(config["output_location"] + "segmentation/{sample}/{window}.txt", sample=SAMPLES, window=[100000]),
        # expand(config["output_location"] + "snv_calls/{sample}/merged.bam", sample=SAMPLES)
        # expand(config["output_location"] + "snv_genotyping/{sample}/{chrom}.vcf", sample=SAMPLES, window=[100000], chrom=config["chromosomes"]),
        # expand(config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz", sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e], window=[100000], ),
        # expand(config["output_location"] + "counts-per-cell/{sample}/{cell}/{window}.txt.gz", sample=SAMPLES, cell=[sub_e for e in list(CELL_PER_SAMPLE.values()) for sub_e in e], window=[100000], ),
        # expand(config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/final.txt", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        expand(config["output_location"] + "haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens}.Rdata", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["input_bam_location"] +  "{sample}/{folder}/{bam}.{chrom}.txt", 
        #         sample=SAMPLES, 
        #         folder=["all", "selected"], 
        #         bam=final_list, 
        #         chrom=config['chromosomes'])

# FIXME : To solve : cell wildcard (dict type) comparatively to others that are list type


################################################################################
# TMP solution by extracting chrom in BAM files                                #
################################################################################


rule simplify_bam_files:
    input:
        bam = config["input_bam_location_raw"] +  "{sample}/{folder}/{bam}.bam"
    output:
        bam_with_header = config["output_location"] + "lite_bam_with_full_header/" +  "{sample}/{folder}/{bam}.bam",
        bam_without_header = config["output_location"] + "lite_bam_with_lite_header/" +  "{sample}/{folder}/{bam}.bam"

    shell:
        """
        cat \
            <(samtools view -H {input.bam} | grep -P "^@HD|^@RG|^\@SQ\\tSN:chr21|^@PG") \
            <(samtools view {input.bam} chr21) |\
        samtools view -bo {output.bam_without_header};
        cat \
            <(samtools view -H {input.bam}) \
            <(samtools view {input.bam} chr21) |\
        samtools view -bo {output.bam_with_header}
        """




################################################################################
# Read counting                                                                #
################################################################################

# CHECKME : exclude file rule useful ?
rule generate_exclude_file_1:
    input:
        bam = expand(config["output_location"] + "lite_bam_with_lite_header/" +  "{sample}/selected/{bam}.bam", sample = SAMPLES[0], bam = BAM_PER_SAMPLE[SAMPLES[0]][0])
    output:
        config["output_location"] + "log/exclude_file.temp"
    log:
        config["output_location"] + "log/generate_exclude_file_1.log"
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -H {input.bam} | awk "/^@SQ/" > {output} 
        """

rule generate_exclude_file_2:
    output:
        config["output_location"] + "log/exclude_file"
    input:
        config["output_location"] + "log/exclude_file.temp"
    params:
        chroms = config["chromosomes"]
    run:
        with open(input[0]) as f:
            with open(output[0],"w") as out:
                for line in f:
                    contig = line.strip().split()[1]
                    contig = contig[3:]
                    # if contig not in params.chroms:
                        # print(contig, file = out)
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
        bam = lambda wc: expand(config["output_location"] + "lite_bam_with_lite_header/" + wc.sample +  "/selected/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        bai = lambda wc: expand(config["output_location"] + "lite_bam_with_lite_header/" + wc.sample +  "/selected/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",

        # excl = config["output_location"] + "log/exclude_file"
    output:
        counts = config["output_location"] + "counts/{sample}/{window}.txt.fixme.gz",
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

rule tmp_filter_mosaic_count_by_chr:
    input:
        config["output_location"] + "counts/{sample}/{window}.txt.fixme.gz"
    output:
        config["output_location"] + "counts/{sample}/{window}.txt.gz"
    run:
        df = pd.read_csv(input[0], compression='gzip', sep='\t')
        df = df.loc[df['chrom'] == 'chr21']
        df.to_csv(output[0], compression='gzip', sep='\t', index=False)



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

# FIXME: no difference observed before/after awk command


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
# TODO : replace R script by integrating directly pandas in the pipeline
# CHECKME : used ???
# DOCME : check if doc is correct
rule prepare_segments:
    """
    rule fct: selection of appropriate segmentation according breakpoint density (k) selected by the user : many : 60%, medium : 40%, few : 20%
    input: mosaic segment output segmentation file
    output: lite file with appropriate k according the quartile defined by the user
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

# TODO : replace awk external file command with something else
rule extract_single_cell_counts:
    """
    rule fct: extract from count the rows coming from the given cell
    input: mosaic count output file for the sample according a given window 
    output: count per cell file for the sample according a given window
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
    rule fct: Same as `rule segmentation` : mosaic segment function but for individual cell
    input: mosaic count splitted by cell produced by `rule extract_single_cell_counts`
    output: Segmentation file for an individual cell
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
    # URGENT : If one bad cell is detected => pipeline stop => need to fix this 

################################################################################
# StrandPhaseR things                                                          #
################################################################################

# DOCME : how to handle when multiple chrom orientation 
"""
RPE1-WT	RPE1WTPE20492	chr10	0	        27300000	WW
RPE1-WT	RPE1WTPE20492	chr10	27300000	110600000	WC
RPE1-WT	RPE1WTPE20492	chr10	110600000	127100000	CC
RPE1-WT	RPE1WTPE20492	chr10	127100000	133797422	WC
"""
"selected_j0.1_s0.5_scedist20"
"""
PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
./utils/detect_strand_states.py \
    --sce_min_distance 500 000 \
    --sce_add_cutoff 20 000 000 \
    --min_diff_jointseg 0.1 \
    --min_diff_singleseg 0.5 \
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
rule segmentation_selection:
    """
    rule fct:
    input: mosaic read counts (txt.gz) & stats info (.info) + joint & sc segmentation 
    output: initial_strand_state used for the following by strandphaser
    """
    input:
        counts=config["output_location"] + "counts/{sample}/{window}.txt.gz",
        jointseg=config["output_location"] + "segmentation/{sample}/{window}.txt",
        singleseg=lambda wc: [config["output_location"] + "segmentation-per-cell/{}/{}/{}.txt".format(wc.sample, cell, wc.window) for cell in CELL_PER_SAMPLE[wc.sample]],
        info=config["output_location"] + "counts/{sample}/{window}.info",
    output:
        jointseg=config["output_location"] + "segmentation2/{sample}/{window,[0-9]+}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        singleseg=config["output_location"] + "segmentation-singlecell/{sample}/{window,[0-9]+}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        strand_states=config["output_location"] + "strand_states/{sample}/{window,[0-9]+}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}/initial_strand_state",
    log:
        config["output_location"] + "log/segmentation_selection/{sample}/{window}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.log"
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
            {input.singleseg} > {log} 2>&1
        """




# TODO : replace R script by integrating directly pandas in the pipeline / potentialy use piped output to following rule ?

rule convert_strandphaser_input:
    """
    rule fct: extract only segmentation with WC orientation 
    input: initial_strand_state file coming from rule segmentation_selection & info file from mosaic count output
    output: filtered TSV file with start/end coordinates of WC-orientated segment to be used by strandphaser
    """
    input:
        states = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/initial_strand_state",
        # URGENT : hard coded 500000 file name ???
        # info   = config["output_location"] + "counts/{sample}/500000.info"
        # FIXME : quick workaround with {window} wc
        info   = config["output_location"] + "counts/{sample}/{window}.info"
    output:
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/strandphaser_input.txt"
    log:
        config["output_location"] + "log/convert_strandphaser_input/{sample}/{window}.{bpdens}.log"
    script:
        "utils/helper.convert_strandphaser_input.R"


# TODO : make something similar to mosaic with C++ dep
# CHECKME : check if possible to write something more snakemak"ic" & compliant with conda/singularity running env
# WARNING : I/O path definition
# WARNING : Try to find a solution to install stranphaser in a conda environment => contact david porubsky to move on the bioconductor ?

# rule install_StrandPhaseR:
#     output:
#         "utils/R-packages/StrandPhaseR/R/StrandPhaseR"
#     log:
#         "log/install_StrandPhaseR.log"
#     shell:
#         """
#         TAR=$(which tar) Rscript utils/install_strandphaser.R > {log} 2>&1
#         """

# TODO : replace by clean config file if possible or by temporary removed file 
rule prepare_strandphaser_config_per_chrom:
    """
    rule fct: prepare config file used by strandphaser
    input: input used only for wildcards : sample, window & bpdens
    output: config file used by strandphaser
    """
    input:
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/initial_strand_state"
    output:
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR.{chrom}.config"
    run:
        with open(output[0], "w") as f:
            print("[General]",                    file = f)
            print("numCPU           = 1",         file = f)
            print("chromosomes      = '" + wildcards.chrom + "'", file = f)
            if (config["paired_end"]):
                print("pairedEndReads   = TRUE",  file = f)
            else:
                print("pairedEndReads   = FALSE", file = f)
            print("min.mapq         = 10",        file = f)
            print("",                             file = f)
            print("[StrandPhaseR]",               file = f)
            print("positions        = NULL",      file = f)
            print("WCregions        = NULL",      file = f)
            print("min.baseq        = 20",       file = f)
            print("num.iterations   = 2",        file = f)
            print("translateBases   = TRUE",     file = f)
            print("fillMissAllele   = NULL",     file = f)
            print("splitPhasedReads = TRUE",     file = f)
            print("compareSingleCells = TRUE",     file = f)
            print("callBreaks       = FALSE",    file = f)
            print("exportVCF        = '", wildcards.sample, "'", sep = "", file = f)
            print("bsGenome         = '", config["R_reference"], "'", sep = "", file = f)



# # TODO : TMP solution 
# # CHECKME : need to check with people if SNP genotyping file is mandatory => will simplify things
# def locate_snv_vcf(wildcards):
#     if "snv_calls" not in config or wildcards.sample not in config["snv_calls"] or config["snv_calls"][wildcards.sample] == "":
#         if "snv_sites_to_genotype" in config and config["snv_sites_to_genotype"] != "" :
#             if os.path.isfile(config["snv_sites_to_genotype"]):
#                 return "snv_genotyping/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
#             else:
#                 return "snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
#         else:
#             return "snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
#     else:
#         return "external_snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)


rule run_strandphaser_per_chrom:
    """
    rule fct: run strandphaser for each chromosome 
    input: strandphaser_input.txt from rule convert_strandphaser_input ; genotyped snv for each chrom by freebayes ; configfile created by rule prepare_strandphaser_config_per_chrom ; bam folder
    output:
    """
    input:
        wcregions    = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/strandphaser_input.txt",
        snppositions = config["output_location"] + "snv_genotyping/{sample}/{chrom}.vcf",
        configfile   = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/StrandPhaseR.{chrom}.config",
        # DOCME : used as an input to call the installation
        # strandphaser = "utils/R-packages/StrandPhaseR/R/StrandPhaseR",
        # strandphaser = config["strandphaser"],
        bamfolder    = config["output_location"] + "lite_bam_with_full_header/" + "{sample}/selected"
    output:
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf"
    log:
        "log/run_strandphaser_per_chrom/{sample}/{window}.{bpdens}/{chrom}.log"
    shell:
        """
        {config[Rscript]} utils/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                {config[output_location]}strand_states/{wildcards.sample}/{wildcards.window}.{wildcards.bpdens}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \

        """

rule compress_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{file}.vcf",
    output:
        vcf="{file}.vcf.gz",
    # log:
    #     "log/compress_vcf/{file}.log"
    shell:
        "bgzip {input.vcf}"


rule index_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{file}.vcf.gz",
    output:
        tbi="{file}.vcf.gz.tbi",
    shell:
        "tabix -p vcf {input.vcf}"

rule merge_strandphaser_vcfs:
    input:
        vcfs=expand(config["output_location"] + "strand_states/{{sample}}/{{window}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"]),
        tbis=expand(config["output_location"] + "strand_states/{{sample}}/{{window}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"]),
    output:
        vcf=config["output_location"] + "phased-snvs/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.vcf.gz"
    log:
        "log/merge_strandphaser_vcfs/{sample}/{window}.{bpdens}.log"
    shell:
        "(bcftools concat -a {input.vcfs} | bcftools view -o {output.vcf} -O z --genotype het --types snps - ) > {log} 2>&1"



rule combine_strandphaser_output:
    input:
        expand(config["output_location"] + "strand_states/{{sample}}/{{window}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt", chrom = config["chromosomes"])
    output:
        config["output_location"] +  "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/strandphaser_output.txt"
    log:
        "log/combine_strandphaser_output/{sample}/{window}.{bpdens}.log"
    shell:
        """
        set +o pipefail
        cat {input} | head -n1 > {output};
        tail -q -n+2 {input} >> {output};
        """


rule convert_strandphaser_output:
    input:
        phased_states  = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/strandphaser_output.txt",
        initial_states = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/initial_strand_state",
        # info           = config["output_location"] + "counts/{sample}/500000_fixed.info"
        info           = config["output_location"] + "counts/{sample}/{window}.info"
    output:
        config["output_location"] + "strand_states/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/final.txt"
    log:
        "log/convert_strandphaser_output/{sample}/{window}.{bpdens}.log"
    script:
        "utils/helper.convert_strandphaser_output.R"



################################################################################
# Haplotagging                                                                 #
################################################################################

rule haplotag_bams:
    input:
        vcf = config["output_location"] + "phased-snvs/{sample}/{window}.{bpdens}.vcf.gz",
        tbi = config["output_location"] + "phased-snvs/{sample}/{window}.{bpdens}.vcf.gz.tbi",
        # WARNING : careful with path
        bam = config["input_bam_location_raw"] + "{sample}/selected/{bam}.bam",
        bai = config["input_bam_location_raw"] + "{sample}/selected/{bam}.bam.bai"
    output:
        bam = config["output_location"] + "haplotag/bam/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/{bam}.bam",
    log:
        config["output_location"] + "log/haplotag_bams/{sample}/{window}.{bpdens}/{bam}.log"
    params:
        ref = config["reference"]
    shell:
        "whatshap haplotag -o {output.bam} -r {params.ref} {input.vcf} {input.bam} > {log} 2>{log}"

rule create_haplotag_segment_bed:
    input:
        segments = config["output_location"] + "segmentation2/{sample}/{size}{what}.{bpdens}.txt",
    output:
        bed = config["output_location"] + "haplotag/bed/{sample}/{size,[0-9]+}{what}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.bed",
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v s={wildcards.size} -f utils/command3.awk {input.segments} > {output.bed}
        """

rule create_haplotag_table:
    input:
        bam = config["output_location"] + "haplotag/bam/{sample}/{window}.{bpdens}/{cell}.bam",
        bai = config["output_location"] + "haplotag/bam/{sample}/{window}.{bpdens}/{cell}.bam.bai",
        bed = config["output_location"] + "haplotag/bed/{sample}/{window}.{bpdens}.bed"
    output:
        tsv = config["output_location"] + "haplotag/table/{sample}/by-cell/haplotag-counts.{cell}.{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.tsv"
    log:
        config["output_location"] + "log/create_haplotag_table/{sample}.{cell}.{window}.{bpdens}.log"
    script:
        "utils/haplotagTable.snakemake.R"

rule merge_haplotag_tables:
    input:
        tsvs = lambda wc: [config["output_location"] + "haplotag/table/{}/by-cell/haplotag-counts.{}.{}.{}.tsv".format(wc.sample,cell,wc.window,wc.bpdens) for cell in BAM_PER_SAMPLE[wc.sample]],
    output:
        tsv = config["output_location"] + "haplotag/table/{sample}/full/haplotag-counts.{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.tsv"
    shell:
        "(head -n1 {input.tsvs[0]} && tail -q -n +2 {input.tsvs}) > {output.tsv}"


rule mosaiClassifier_calc_probs:
    input:
        counts = config["output_location"] + "counts/{sample}/{window}.txt.gz",
        info   = config["output_location"] + "counts/{sample}/{window}.info",
        states = config["output_location"] + "strand_states/{sample}/{window}.{bpdens}/final.txt",
        bp     = config["output_location"] + "segmentation2/{sample}/{window}.{bpdens}.txt"
    output:
        output = config["output_location"] + "sv_probabilities/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/probabilities.Rdata"
    log:
        config["output_location"] + "log/mosaiClassifier_calc_probs/{sample}/{window}.{bpdens}.log"
    script:
        "utils/mosaiClassifier.snakemake.R"

rule create_haplotag_likelihoods:
    input:
        haplotag_table = config["output_location"] + "haplotag/table/{sample}/full/haplotag-counts.{window}.{bpdens}.tsv",
        sv_probs_table = config["output_location"] + "sv_probabilities/{sample}/{window}.{bpdens}/probabilities.Rdata",
    output: 
        config["output_location"] + "haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.Rdata"
    log:
        config["output_location"] + "log/create_haplotag_likelihoods/{sample}.{window}.{bpdens}.log"
    params:
        chroms = config["chromosomes"]
    shell:
        "/home/tweber/.conda/envs/strandseqnation/bin/Rscript afac/haplotagProbs.snakemake.R {input.haplotag_table} {input.sv_probs_table} {output}"




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
        lambda wc: expand(config["output_location"] + "lite_bam_with_full_header/" + wc.sample + "/all/{bam}.bam", bam = ALLBAMS_PER_SAMPLE[wc.sample]) if wc.sample in ALLBAMS_PER_SAMPLE else "FOOBAR",
    output:
        config["output_location"] + "snv_calls/{sample}/merged.unsorted.bam"
    log:
        config["output_location"] + "log/mergeBams/{sample}.log"
    threads:
        4
    shell:
        # FIXME : Samtools 1.10 from Conda env not working ; 1.9 from Seneca working > change it into conda env yml file
        config["samtools"] + " merge -@ {threads} {output} {input} 2>&1 > {log}"
        # "samtools" + " merge -@ {threads} {output} {input} 2>&1 > {log}"

rule sort_bam:
    input:
        config["output_location"] + "snv_calls/{sample}/merged.unsorted.bam"
    output:
        config["output_location"] + "snv_calls/{sample}/merged.bam"
    shell:
        config["samtools"] + " sort {input} -o {output}"


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
            --include "QUAL>=10" \
        > {output.vcf}) 2> {log}
        """
