import math
from collections import defaultdict

configfile: "Snake.config.yaml"

import os, sys

# print(os.listdir(os.getcwd()))
# print(os.listdir("bam"))

SAMPLE,BAM = glob_wildcards(config['input_bam_location'] + "{sample}/selected/{bam}.bam")

print(SAMPLE)
print(BAM)

SAMPLES = sorted(set(SAMPLE))

print(SAMPLES)
from pprint import pprint


CELL_PER_SAMPLE= defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)


for sample,bam in zip(SAMPLE,BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace('.sort.mdup',''))

pprint(BAM_PER_SAMPLE)
pprint(CELL_PER_SAMPLE)


ALLBAMS_PER_SAMPLE = defaultdict(list)
for sample in SAMPLES:
    ALLBAMS_PER_SAMPLE[sample] = glob_wildcards(config['input_bam_location'] + "{}/all/{{bam}}.bam".format(sample)).bam

pprint(ALLBAMS_PER_SAMPLE)



print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print("  {}:\t{} cells\t {} selected cells".format(s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])))


exit()

import os.path

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
        expand("counts/{sample}/{window}_fixed.txt.gz", sample=SAMPLES, window=[100000])


################################################################################
# Read counting                                                                #
################################################################################

rule generate_exclude_file_1:
    output:
        temp("log/exclude_file.temp")
    input:
        bam = expand("bam/{sample}/selected/{bam}.bam", sample = SAMPLES[0], bam = BAM_PER_SAMPLE[SAMPLES[0]][0])
    log:
        "log/generate_exclude_file_1.log"
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -H {input.bam} | awk "/^@SQ/" > {output} 
        """

rule generate_exclude_file_2:
    output:
        "log/exclude_file"
    input:
        "log/exclude_file.temp"
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


rule mosaic_count_fixed:
    input:
        bam = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        bai = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        excl = "log/exclude_file"
    output:
        counts = "counts/{sample}/{window}_fixed.txt.gz",
        info   = "counts/{sample}/{window}_fixed.info"
    log:
        "log/{sample}/mosaic_count_fixed.{window}.log"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        echo mosaic_count_fixed && 
        {params.mc_command} count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {wildcards.window} \
            {input.bam} 
        """

# rule extract_single_cell_counts:
#     input:
#         "counts/{sample}/{window}_{file_name}.txt.gz"
#     output:
#         "counts-per-cell/{sample}/{cell}/{window,[0-9]+}_{file_name}.txt.gz"
#     shell:
#         """
#         # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
#         zcat {input} | awk -v name={wildcards.cell} -f utils/command1.awk | gzip > {output}
#         """


# ################################################################################
# # Plots                                                                        #
# ################################################################################

# rule plot_mosaic_counts:
#     input:
#         counts = "counts/{sample}/{file_name}.txt.gz",
#         info   = "counts/{sample}/{file_name}.info"
#     output:
#         "plots/{sample}/{file_name}.pdf"
#     log:
#         "log/plot_mosaic_counts/{sample}/{file_name}.log"
#     params:
#         plot_command = "Rscript " + config["plot_script"]
#     shell:
#         """
#         {params.plot_command} {input.counts} {input.info} {outp
#         """

