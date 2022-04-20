# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher

import math
from collections import defaultdict

configfile: "config/Snake.config_embl.yaml"
import pandas as pd
import os, sys
from pprint import pprint
import pysam
from tqdm import tqdm

# TODO I/O : Function to define inputs ; simplify list/dict system
# TODO Use remote file system to download example files


def handle_input_data(thisdir, exclude_list=list):
    """
        
    """
    # Parsing folder and retrieve only files with .bam extension
    data = [(r,file.replace('.bam', '')) for r, d, f in os.walk(thisdir) for file in f if ".bam" in file and ".bai" not in file]
    
    # Building pandas df based on folder structure
    df = pd.DataFrame(data,columns=['Folder','File'])

    # Defining cols
    df['all/selected'] = df['Folder'].apply(lambda r: r.split('/')[-1])
    df['Sample'] = df['Folder'].apply(lambda r: r.split('/')[-2])
    df['Cell'] = df['File'].apply(lambda r: r.split('.')[0])
    df['Full_path'] = df['Folder'] + "/" + df['File'] + ".bam"

    # Filtering based on exclude list defined
    df_config_files = df.loc[~df['Cell'].isin(exclude_list)]

    # Export dicts
    SAMPLES = sorted(df_config_files.Sample.unique().tolist())
    BAM_PER_SAMPLE = df_config_files.loc[df_config_files['all/selected'] == "selected"].groupby('Sample')['File'].apply(list).to_dict()
    CELL_PER_SAMPLE = df_config_files.loc[df_config_files['all/selected'] == "selected"].groupby('Sample')['Cell'].apply(list).to_dict()
    ALLBAMS_PER_SAMPLE = df_config_files.loc[df_config_files['all/selected'] == "all"].groupby('Sample')['File'].apply(list).to_dict()

    return SAMPLES, BAM_PER_SAMPLE, CELL_PER_SAMPLE, ALLBAMS_PER_SAMPLE, df_config_files


def check_bam_header(bam_file_path):
    """
        
    """

    # Get BAM file header with pysam
    h = pysam.view("-H", bam_file_path)
    h = [e.split("\t") for e in h.split("\n")]
    sm_tag_list = list(set([sub_e.replace("SM:", "") for e in h for sub_e in e if "SM:" in sub_e]))

    # Folder name based on path
    folder_name = bam_file_path.split("/")[-3]

    # Assertions
    assert len(sm_tag_list) == 1, "Two different SM tags in the header of BAM file {}".format(bam_file_path)
    assert sm_tag_list[0] == folder_name, 'Folder name "{}" must correspond to SM tag in BAM file "{}"'.format(folder_name, bam_file_path)


# FIXME : tmp solution to remove bad cells => need to fix this with combination of ASHLEYS ?
# TODO : other solution by giving in config file, CLI input ?

exclude_list = ['BM510x3PE20490']

SAMPLES, BAM_PER_SAMPLE, CELL_PER_SAMPLE, ALLBAMS_PER_SAMPLE, df_config_files = handle_input_data(thisdir=config["input_bam_location"], exclude_list=exclude_list)

print(df_config_files)
print(df_config_files['Full_path'][0])

tqdm.pandas(desc="Checking if BAM SM tags correspond to folder names")
df_config_files["Full_path"].progress_apply(check_bam_header, )

print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print("  {}:\t{} cells\t {} selected cells".format(s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])))



METHODS = [
    "simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE",
]

# # FIXME : move to yaml/json settings or to something else
BPDENS = [
    "selected_j{}_s{}_scedist{}".format(joint, single, scedist) for joint in [0.1] for single in [0.5] for scedist in [20]
]


rule all:
    input:
        # expand(config["output_location"] + "counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
        expand(config["output_location"] + "counts/{sample}.txt.gz", sample=SAMPLES),





# FIXME : To solve : cell wildcard (dict type) comparatively to others that are list type




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
        bam = df_config_files['Full_path'][0]
    output:
        config["output_location"] + "log/exclude_file"
    params:
        chroms = config["chromosomes"]
    run:
        # READ BAM FILE HEADER OF FIRST BAM IN THE PANDAS DF
        h = pysam.view("-H", input[0])
        h = [e.split("\t") for e in h.split("\n") if "@SQ" in e]

        # CONVERT TO PANDAS DF
        df_h = pd.DataFrame(h, columns=["TAG", "Contig", "LN"])

        # PROCESS CONTIGS
        output_h = pd.DataFrame(df_h["Contig"].str.replace("SN:", ""))
        output_h = output_h.loc[~output_h["Contig"].isin(params.chroms)]
        
        # EXPORT
        output_h["Contig"].to_csv(output[0], index=False, sep='\t', header=False)

                        
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
        bai = lambda wc: expand(config["input_bam_location"] + wc.sample +  "/selected/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        excl = config["output_location"] + "log/exclude_file"
    output:
        counts = config["output_location"] + "counts/{sample}.txt.gz",
        info   = config["output_location"] + "counts/{sample}.info"
    log:
        config["output_location"] + "log/{sample}/mosaic_count.log"
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
            {input.bam} 
        > {log} 2>&1
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


################################################################################
# REGENOTYPE SNV                                                               #
################################################################################

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


################################################################################
# StrandPhaseR things                                                          #
################################################################################


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
        bamfolder    = config["input_bam_location"] + "{sample}/selected"
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
        bam = config["input_bam_location"] + "{sample}/selected/{bam}.bam",
        bai = config["input_bam_location"] + "{sample}/selected/{bam}.bam.bai"
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





################################################################################
# MosaiClassifier                                                              #
################################################################################

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
        haplotag_table = config["output_location"] + 'haplotag/table/{sample}/full/haplotag-counts.{window}.{bpdens}.tsv',
        sv_probs_table = config["output_location"] + 'sv_probabilities/{sample}/{window}.{bpdens}/probabilities.Rdata',
    output: 
        config["output_location"] + 'haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.Rdata'
    log:
        config["output_location"] + "log/create_haplotag_likelihoods/{sample}.{window}.{bpdens}.log"
    script:
        "utils/haplotagProbs.snakemake.R"

rule mosaiClassifier_make_call:
    input:
        probs = config["output_location"] + 'haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens}.Rdata'
    output:
        config["output_location"] + "sv_calls/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filterFALSE.txt"
    params:
        minFrac_used_bins = 0.8
    log:
        config["output_location"] + "log/mosaiClassifier_make_call/{sample}/{window}.{bpdens}.llr{llr}.poppriors{pop_priors}.haplotags{use_haplotags}.gtcutoff{gtcutoff}.regfactor{regfactor}.log"
    script:
        "utils/mosaiClassifier_call.snakemake.R"



# CHECKME : check if still useful ?
rule mosaiClassifier_make_call_biallelic:
    input:
        probs = config["output_location"] + "sv_probabilities/{sample}/{window}.{bpdens}/probabilities.Rdata"
    output:
        config["output_location"] + "sv_calls/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/biAllelic_llr{llr}.txt"
    log:
        config["output_location"] + "log/mosaiClassifier_make_call_biallelic/{sample}/{window}.{bpdens}.{llr}.log"
    script:
        "utils/mosaiClassifier_call_biallelic.snakemake.R"



################################################################################
# PostProcessing                                                               #
################################################################################


# DOCME : perl in conda

rule postprocessing_filter:
    input: 
        calls = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.txt"
    output: 
        calls = config["output_location"] + "postprocessing/filter/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.txt"
    shell:
        'utils/filter_MosaiCatcher_calls.pl {input.calls}  > {output.calls}'

rule postprocessing_merge:
    input: 
        calls = config["output_location"] + "postprocessing/filter/{sample}/{window}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.txt"
    output: 
        calls = config["output_location"] + "postprocessing/merge/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.txt"
    shell:
        'utils/group_nearby_calls_of_same_AF_and_generate_output_table.pl {input.calls}  > {output.calls}'


rule postprocessing_sv_group_table:
    input: 
        calls = config["output_location"] + "postprocessing/merge/{sample}/{window}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.txt"
    output: 
        grouptrack = config["output_location"] + "postprocessing/group-table/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.tsv"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        utils/create-sv-group-track.py {input.calls}  > {output.grouptrack}
        """



rule filter_calls:
    input: 
        inputcalls = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.txt",
        mergedcalls = config["output_location"] + "postprocessing/merge/{sample}/{window}.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.txt",
    output: 
        calls = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filterTRUE.txt"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        utils/apply_filter.py {input.inputcalls} {input.mergedcalls} > {output.calls}
        """



rule call_complex_regions:
    input:
        calls  = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/{method}_filter{filter}.txt",
    output:
        complex = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/{method}_filter{filter}.complex.tsv",
    log:
        config["output_location"] + "log/call_complex_regions/{sample}/{window}.{bpdens}.{method}_filter{filter}.log"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        utils/call-complex-regions.py \
        --merge_distance 5000000 \
        --ignore_haplotypes \
        --min_cell_count 2 {input.calls} > {output.complex} 2>{log}
        """



################################################################################
# Summary statistics on sv calls                                               #
################################################################################


rule summary_statistics:
    input:
        segmentation = config["output_location"] + 'segmentation2/{sample}/{window}.{bpdens}.txt',
        strandstates = config["output_location"] + 'strand_states/{sample}/{window}.{bpdens}/initial_strand_state',
        sv_calls = config["output_location"] + 'sv_calls/{sample}/{window}.{bpdens}/{method}_filter{filter}.txt',
        complex = config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/{method}_filter{filter}.complex.tsv",
        merged = config["output_location"] + "postprocessing/merge/{sample}/{window}.{bpdens}/{method}.txt",
    output:
        tsv = config["output_location"] + 'stats/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/{method}_filter{filter,(TRUE|FALSE)}.tsv',
    log:
        config["output_location"] + 'log/summary_statistics/{sample}/{window}.{bpdens}/{method}_filter{filter}.log'
    run:
        p = []
        try:
            f = config["ground_truth_clonal"][wildcards.sample]
            if len(f) > 0:
                p.append('--true-events-clonal')
                p.append(f)
        except KeyError:
            pass
        try:
            f = config["ground_truth_single_cell"][wildcards.sample]
            if len(f) > 0:
                p.append('--true-events-single-cell')
                p.append(f)
        except KeyError:
            pass
        if wildcards.filter == 'TRUE':
            p.append('--merged-file')
            p.append(input.merged)
        additional_params = ' '.join(p)
        shell('utils/callset_summary_stats.py --segmentation {input.segmentation} --strandstates {input.strandstates} --complex-regions {input.complex} {additional_params} {input.sv_calls}  > {output.tsv} ')

rule aggregate_summary_statistics:
    input:
        tsv=expand(config["output_location"] + "stats/{{sample}}/{window}.{bpdens}/{method}.tsv", window = [100000], bpdens = BPDENS, method = METHODS),
    output:
        tsv=config["output_location"] + "stats-merged/{sample}/stats.tsv"
    shell:
        "(head -n1 {input.tsv[0]} && (tail -n1 -q {input.tsv} | sort -k1) ) > {output}"
    

# CHECKME : to check & see if it's working
################################################################################
# Ploidy estimation                                                            #
################################################################################
# TODO : merge into one file by sample 
rule estimate_ploidy:
    input:
        config["output_location"] + "counts/{sample}/100000.txt.gz"
    output:
        config["output_location"] + "ploidy/{sample}/ploidy.{chrom}.txt"
    log:
        config["output_location"] + "log/estimate_ploidy/{sample}/{chrom}.log"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        python utils/ploidy-estimator.py --chromosome {wildcards.chrom} {input} > {output} 
        """


################################################################################
# Plots                                                                        #
################################################################################
# FIXME : Missing plots in final PDF ; R script + inputs to check
# CHECKME : check if possible to switch from PDF to svg (or both) to produce lighter files
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
        plot_command = config["Rscript"] + " " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output} > {log} 2>&1
        """


rule generate_halo_json:
    input:
        counts = config["output_location"] + "counts/{sample}/{windows}.txt.gz",
    output:
        json = config["output_location"] + "halo/{sample}/{windows}.json.gz",
    log:
        config["output_location"] + "log/generate_halo_json/{sample}/{windows}.{windows}.log"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        (./utils/counts_to_json.py {input.counts} | gzip > {output.json}) 
        """


rule plot_SV_calls:
    input:
        counts = config["output_location"] + "counts/{sample}/{windows}.txt.gz",
        calls  = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.txt",
        complex = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.complex.tsv",
        strand = config["output_location"] + "strand_states/{sample}/{windows}.{bpdens}/final.txt",
        segments = config["output_location"] + "segmentation2/{sample}/{windows}.{bpdens}.txt",
        scsegments = config["output_location"] + "segmentation-singlecell/{sample}/{windows}.{bpdens}.txt",
        grouptrack = config["output_location"] + "postprocessing/group-table/{sample}/{windows}.{bpdens}/{method}.tsv",
    output:
        config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_calls/{method}_filter{filter,(TRUE|FALSE)}.{chrom}.pdf"
    log:
        config["output_location"] + "log/plot_SV_calls/{sample}/{windows}.{bpdens}.{method}_filter{filter}.{chrom}.log"
    shell:
        """
        {config[Rscript]} utils/plot-sv-calls.R \
            segments={input.segments} \
            singlecellsegments={input.scsegments} \
            strand={input.strand} \
            complex={input.complex} \
            groups={input.grouptrack} \
            calls={input.calls} \
            {input.counts} \
            {wildcards.chrom} \
            {output} > {log} 2>&1
        """



rule plot_SV_consistency_barplot:
    input:
        sv_calls  = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt",
    output:
        barplot_bypos = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_consistency/{method}.consistency-barplot-bypos.pdf",
        barplot_byaf = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_consistency/{method}.consistency-barplot-byaf.pdf",
    log:
        config["output_location"] + "log/plot_SV_consistency/{sample}/{windows}.{bpdens}.{method}.log"
    shell:
        """
        {config[Rscript]} utils/sv_consistency_barplot.snakemake.R
        """




rule plot_clustering:
    input:
        sv_calls  = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt",
        binbed = "utils/bin_200kb_all.bed",
    output:
        position = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_clustering/{method}-position.pdf",
        chromosome = config["output_location"] + "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_clustering/{method}-chromosome.pdf",
    log:
        config["output_location"] + "log/plot_clustering/{sample}/{windows}.{bpdens}.{method}.log"
    shell:
        """
        {config[Rscript]} utils/plot-clustering.snakemake.R {input.sv_calls} {input.binbed} {output.position} {output.chromosome}
        """



################################################################################
# UTILS                                                                        #
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



"""
# RXIV

rule all:
    input:
        expand(config["output_location"] + "plots/{sample}/{window}.pdf", sample = SAMPLES, window = [100000]),

        # expand(config["output_location"] + "counts/{sample}/{window}.txt.gz", sample=SAMPLES, window=[100000]),
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
        # expand(config["output_location"] + "haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens}.Rdata", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["output_location"] + "sv_probabilities/{sample}/{window}.{bpdens}/probabilities.Rdata", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/biAllelic_llr4.txt", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/biAllelic_llr4.complex.tsv", sample=SAMPLES, window=[100000], bpdens=BPDENS, chrom=config["chromosomes"]),
        # expand(config["output_location"] + "postprocessing/merge/{sample}/{window}.{bpdens}/{method}.txt",
        #        sample = SAMPLES,
        #        window = [100000],
        #        bpdens = BPDENS,
        #        method = list(set(m.replace('_filterTRUE','').replace('_filterFALSE','') for m in METHODS))),
        # expand(config["output_location"] + "stats-merged/{sample}/stats.tsv", sample = SAMPLES),

        # expand(config["input_bam_location"] +  "{sample}/{folder}/{bam}.{chrom}.txt", 
        #         sample=SAMPLES, 
        #         folder=["all", "selected"], 
        #         bam=final_list, 
        #         chrom=config['chromosomes'])
        expand(config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/plots/sv_calls/{method}.{chrom}.pdf",
               sample = SAMPLES,
               chrom = config["chromosomes"],
               window = [100000],
               bpdens = BPDENS,
               method = METHODS),
        # expand("ploidy/{sample}/ploidy.{chrom}.txt", sample = SAMPLES, chrom = config["chromosomes"]),
        expand(config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/plots/sv_consistency/{method}.consistency-barplot-{plottype}.pdf",
               sample = SAMPLES,
               window = [100000],
               bpdens = BPDENS,
               method = METHODS,
               plottype = ["byaf","bypos"]),
        expand(config["output_location"] + "sv_calls/{sample}/{window}.{bpdens}/plots/sv_clustering/{method}-{plottype}.pdf",
               sample = SAMPLES,
               window = [100000],
               bpdens = BPDENS,
               method = METHODS,
               plottype = ["position","chromosome"]),
        expand(config["output_location"] + "halo/{sample}/{window}.json.gz",
               sample = SAMPLES,
               window = [100000]),
        expand(config["output_location"] + "ploidy/{sample}/ploidy.{chrom}.txt", sample = SAMPLES, chrom = config["chromosomes"]),

        # expand("stats-merged/{sample}/stats.tsv", sample = SAMPLES),
        # expand("postprocessing/merge/{sample}/{window}.{bpdens}/{method}.txt",
        #        sample = SAMPLES,
        #        window = [100000],
        #        bpdens = BPDENS,
        #        method = list(set(m.replace('_filterTRUE','').replace('_filterFALSE','') for m in METHODS))),


"""