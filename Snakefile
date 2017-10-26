configfile: "Snake.config.json"

BAM, = glob_wildcards("bam/{bam}.bam")

# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher

rule all:
    input:
        expand("plots/" + config["sample"] + ".{window}_fixed.pdf", window = [50000, 100000, 200000, 500000]),
        expand("plots/" + config["sample"] + ".{window}_variable.pdf", window = [50000, 100000]),
        expand("segmentation2/" + config["sample"] + ".{window}_fixed.{bpdens}.txt", window = [50000, 100000, 200000, 500000], bpdens = ["few","many"]),
        expand("segmentation2/" + config["sample"] + ".{window}_variable.{bpdens}.txt", window = [50000, 100000], bpdens = ["few","many"]),
        "strand_states/" + config["sample"] + ".txt"

rule plot_mosaic_counts:
    input:
        counts = "counts/" + config["sample"] + ".{file_name}.txt.gz",
        info   = "counts/" + config["sample"] + ".{file_name}.info"
    output:
        "plots/" + config["sample"] + ".{file_name}.pdf"
    params: 
        plot_command = "Rscript " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output}
        """

rule mosaic_count_fixed:
    input:
        bam = expand("bam/{bam}.bam", bam = BAM),
        bai = expand("bam/{bam}.bam.bai", bam = BAM)
    output:
        counts = "counts/" + config["sample"] + "_{window}_fixed.txt.gz",
        info   = "counts/" + config["sample"] + "_{window}_fixed.info"
    params:
        mc_command = config["mosaicatcher"],
        mc_exclfile = config["exclude_file"]
    shell:
        """
        {params.mc_command} count \
            -o {output.counts} \
            -i {output.info} \
            -x {params.mc_exclfile} \
            -w {wildcards.window} \
            {input.bam}
        """


rule mosaic_count_variable:
    input:
        bam = expand("bam/{bam}.bam", bam = BAM),
        bai = expand("bam/{bam}.bam.bai", bam = BAM),
        bed = lambda wc: config["variable_bins"][str(wc.window)]
    output:
        counts = "counts/" + config["sample"] + "_{window}_variable.txt.gz",
        info   = "counts/" + config["sample"] + "_{window}_variable.info"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} count \
            -o {output.counts} \
            -i {output.info} \
            -b {input.bed} \
            {input.bam}
        """




rule determine_strand_states:
    input:
        "counts/" + config["sample"] + "_500000_fixed.txt.gz"
    output:
        "strand_states/" + config["sample"] + ".txt"
    params:
        sce_command = "Rscript " + config["sce_script"]
    shell:
        """
        {params.sce_command} {input} {output}
        """




################################################################################
# Segmentation                                                                 #
################################################################################

rule segmentation:
    input:
        "counts/" + config["sample"] + ".{file_name}.txt.gz"
    output:
        "segmentation/" + config["sample"] + ".{file_name}.txt"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} segment \
        -o {output} \
        {input}
        """

# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    input:
        "segmentation/" + config["sample"] + ".{windows}.txt"
    output:
        "segmentation2/" + config["sample"] + ".{windows}.{bpdens}.txt"
    params:
        quantile = lambda wc: config["bp_density"][wc.bpdens]
    script:
        "utils/helper.prepare_segments.R"



# Run SV classification
rule run_sv_classification:
    input:
        counts = "counts/" + config["sample"] + ".{windows}.txt.gz",
        info   = "counts/" + config["sample"] + ".{windows}.info",
        states = "strand_states/" + config["sample"] + ".txt",
        bp     = "segmentation2/" + config["sample"] + ".{windows}.{bpdens}.txt"
    output:
        outdir = "sv_probabilities/{file_name}.{bpdens}/",
        outfile1 = "sv_probabilities/{file_name}.{bpdens}/YYY"
    params:
        class_dir     = config["class_dir"],
        class_command = "Rscript " + config["class_script"]
    shell:
        """
        cd {params.class_dir}
        {params.class_command} \
            binRCfile={input.counts} \
            BRfile={input.bp} \
            infoFile={input.info} \
            stateFile={input.states} \
            K=22 \
            maxCN= 4 \
            outputDir={output.outdir}
        """



