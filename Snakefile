configfile: "Snake.config.json"

SAMPLE,BAM = glob_wildcards("bam/{sample}/{bam}.bam")
BAM_PER_SAMPLE = dict([(s,[]) for s in SAMPLE])
for i in range(len(SAMPLE)):
    BAM_PER_SAMPLE[SAMPLE[i]].append(BAM[i])
print("Detected {} samples:".format(len(set(SAMPLE))))
for s in set(SAMPLE):
    print("  {}:\t{} cells".format(s, len(BAM_PER_SAMPLE[s])))


# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher

rule all:
    input:
        expand("plots/{sample}/{window}_fixed.pdf", sample = SAMPLE, window = [50000, 100000, 200000, 500000]),
        expand("plots/{sample}/{window}_variable.pdf", sample = SAMPLE, window = [50000, 100000]),
        expand("segmentation2/{sample}/{window}_fixed.{bpdens}.txt", sample = SAMPLE,
               window = [50000, 100000, 200000, 500000], bpdens = ["few","medium","many"]),
        #expand("segmentation2/{sample}/{window}_variable.{bpdens}.chr1.txt", sample = SAMPLE,
        #       window = [50000, 100000], bpdens = ["few","medium","many"]),
        expand("strand_states/{sample}/final.txt", sample = SAMPLE),
        expand("sv_calls/{sample}/{window}_fixed.{bpdens}.SV_probs.chr1.pdf", sample = SAMPLE,
               window = [50000, 100000, 200000, 500000], bpdens = ["few","medium","many"]),
        #expand("sv_calls/{sample}/{window}_variable.{bpdens}.SV_probs.chr1.pdf", sample = SAMPLE,
        #       window = [50000, 100000], bpdens = ["few","medium","many"])



################################################################################
# Plots                                                                        #
################################################################################

rule plot_mosaic_counts:
    input:
        counts = "counts/{sample}/{file_name}.txt.gz",
        info   = "counts/{sample}/{file_name}.info"
    output:
        "plots/{sample}/{file_name}.pdf"
    log:
        "log/plot_mosaic_counts/{sample}.{file_name}.txt"
    params:
        plot_command = "Rscript " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output} > {log} 2>&1
        """

rule plot_SV_calls:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        probs  = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.txt"
    output:
        dynamic("sv_calls/{sample}/{windows}.{bpdens}.SV_probs.{chrom}.pdf")
    log:
        "log/plot_SV_call/{sample}.{windows}.{bpdens}.txt"
    params:
        plot_command = "Rscript " + config["sv_plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.probs} {output}  > {log} 2>&1
        """


################################################################################
# Read counting                                                                #
################################################################################

rule generate_exclude_file_1:
    output:
        temp("log/exclude_file.temp")
    input:
        bam = expand("bam/{sample}/{bam}.bam", sample = SAMPLE[0], bam = BAM[0]),
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -H {input.bam} | awk '/^@SQ/ {{print substr($2,4)}}' > {output}
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
                    if line.strip() not in params.chroms:
                        print(line.strip(), file = out)



rule mosaic_count_fixed:
    input:
        bam = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]),
        bai = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]),
        excl = "log/exclude_file"
    output:
        counts = "counts/{sample}/{window}_fixed.txt.gz",
        info   = "counts/{sample}/{window}_fixed.info"
    log:
        "log/{sample}/mosaic_count_fixed.{window}.txt"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} count \
            --verbose \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {wildcards.window} \
            {input.bam} \
        > {log} 2>&1
        """


rule mosaic_count_variable:
    input:
        bam = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]),
        bai = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]),
        bed = lambda wc: config["variable_bins"][str(wc.window)],
        excl = "log/exclude_file"
    output:
        counts = "counts/{sample}/{window}_variable.txt.gz",
        info   = "counts/{sample}/{window}_variable.info"
    log:
        "log/{sample}/mosaic_count_variable.{window}.txt"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        echo "NOTE: Exclude file not used in variable-width bins"
        {params.mc_command} count \
            --verbose \
            -o {output.counts} \
            -i {output.info} \
            -b {input.bed} \
            {input.bam} \
        > {log} 2>&1
        """






################################################################################
# Segmentation                                                                 #
################################################################################

rule segmentation:
    input:
        "counts/{sample}/{file_name}.txt.gz"
    output:
        "segmentation/{sample}/{file_name}.txt"
    log:
        "log/{sample}/segmentation.{file_name}.txt"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        {params.mc_command} segment \
        -o {output} \
        {input} > {log} 2>&1
        """

# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    input:
        "segmentation/{sample}/{windows}.txt"
    output:
        "segmentation2/{sample}/{windows}.{bpdens}.txt"
    log:
        "log/{sample}/prepare_segments.{windows}.{bpdens}.txt"
    params:
        quantile = lambda wc: config["bp_density"][wc.bpdens]
    script:
        "utils/helper.prepare_segments.R"


################################################################################
# SV classification                                                            #
################################################################################

rule install_MaRyam:
    output:
        "utils/R-packages2/MaRyam/R/MaRyam"
    log:
        "log/install_MaRyam.log"
    shell:
        """
        Rscript utils/install_maryam.R > {log} 2>&1
        """

rule run_sv_classification:
    input:
        maryam = "utils/R-packages2/MaRyam/R/MaRyam",
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        outdir = "sv_probabilities/{sample}/{windows}.{bpdens}/",
        out1   = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table"
    log:
        "log/{sample}/run_sv_classification.{windows}.{bpdens}.txt"
    params:
        windowsize    = lambda wc: wc.windows.split("_")[0]
    shell:
        """
        set -x
        # set haplotypeInfo if phasing info is available
        Rscript utils/MaRyam_pipeline.R \
                binRCfile={input.counts} \
                BRfile={input.bp} \
                infoFile={input.info} \
                stateFile={input.states} \
                outputDir={output.outdir} \
                bin.size={params.windowsize} \
                K=22 \
                maximumCN=4 \
                utils/R-packages2/ > {log} 2>&1
        """

rule convert_SVprob_output:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table",
        info  = "counts/{sample}/{windows}.info"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.txt"
    params:
        sample_name = "{wildcards.sample}"
    log:
        "log/{sample}/convert_SVprob_output.{windows}.{bpdens}.txt"
    script:
        "utils/helper.convert_svprob_output.R"


################################################################################
# Strand states & phasing                                                      #
################################################################################

rule determine_initial_strand_states:
    input:
        "counts/{sample}/500000_fixed.txt.gz"
    output:
        "strand_states/{sample}/intitial_strand_state"
    log:
        "log/{sample}/determine_initial_strand_states.txt"
    params:
        sce_command = "Rscript " + config["sce_script"]
    shell:
        """
        {params.sce_command} {input} {output}
        """

# Strandphaser needs a different input format which contains the path names to
# the bam files. This rule extracts this information and prepares an input file.
rule convert_strandphaser_input:
    input:
        states = "strand_states/{sample}/intitial_strand_state",
        info   = "counts/{sample}/500000_fixed.info"
    output:
        "strand_states/{sample}/strandphaser_input.txt"
    log:
        "log/{sample}/convert_strandphaser_input.txt"
    script:
        "utils/helper.convert_strandphaser_input.R"

rule install_StrandPhaseR:
    output:
        "utils/R-packages/StrandPhaseR/R/StrandPhaseR"
    log:
        "log/strandphaser-install.log"
    shell:
        """
        Rscript utils/install_strandphaser.R > {log} 2>&1
        """

rule prepare_strandphaser_config_per_chrom:
    input:
        "strand_states/{sample}/intitial_strand_state"
    output:
        "strand_states/{sample}/StrandPhaseR.{chrom}.config"
    run:
        with open(output[0], "w") as f:
            print("[General]",                    file = f)
            print("numCPU           = 1",         file = f)
            print("chromosomes      = '" + wildcards.chrom + "'", file = f)
            print("pairedEndReads   = TRUE",      file = f)
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
            print("exportVCF        = '" + wildcards.sample + ".txt'", sep = "", file = f)
            print("bsGenome         = '", config["R_reference"], "'", sep = "", file = f)


def locate_snv_vcf(wildcards):
    if "snv_calls" not in config or wildcards.sample not in config["snv_calls"] or config["snv_calls"][wildcards.sample] == "":
        return "snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
    else:
        return "external_snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)

rule run_strandphaser_per_chrom:
    input:
        wcregions    = "strand_states/{sample}/strandphaser_input.txt",
        snppositions = locate_snv_vcf,
        configfile   = "strand_states/{sample}/StrandPhaseR.{chrom}.config",
        strandphaser = "utils/R-packages/StrandPhaseR/R/StrandPhaseR",
        bamfolder    = "bam/{sample}/"
    output:
        "strand_states/{sample}/StrandPhaseR_analysis.{chrom}/phased_haps.txt"
    log:
        "log/{sample}/run_strandphaser.{chrom}.txt"
    shell:
        """
        Rscript utils/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                strand_states/{wildcards.sample}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \
                > {log} 2>&1
        """



rule combine_strandphaser_output:
    input:
        expand("strand_states/{{sample}}/StrandPhaseR_analysis.{chrom}/phased_haps.txt",
                chrom = config["chromosomes"])
    output:
        "strand_states/{sample}/strandphaser_output.txt"
    shell:
        """
        set +o pipefail
        cat {input} | head -n1 > {output};
        for x in {input}; do tail -n+2 $x >> {output}; done;
        """


rule convert_strandphaser_output:
    input:
        phased_states  = "strand_states/{sample}/strandphaser_output.txt",
        initial_states = "strand_states/{sample}/intitial_strand_state",
        info           = "counts/{sample}/500000_fixed.info"
    output:
        "strand_states/{sample}/final.txt"
    log:
        "log/{sample}/convert_strandphaser_output.txt"
    script:
        "utils/helper.convert_strandphaser_output.R"



################################################################################
# Call SNVs                                                                    #
################################################################################

rule mergeBams:
    input:
        lambda wc: expand("bam/" + wc.sample + "/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample])
    output:
        "snv_calls/{sample}/merged.bam"
    shell:
        config["samtools"] + " merge {output} {input}"

rule indexMergedBam:
    input:
        "snv_calls/{sample}/merged.bam"
    output:
        "snv_calls/{sample}/merged.bam.bai"
    shell:
        config["samtools"] + " index {input}"


rule call_SNVs_bcftools_chrom:
    input:
        fa    = config["reference"],
        bam   = "snv_calls/{sample}/merged.bam",
        bai   = "snv_calls/{sample}/merged.bam.bai"
    output:
        "snv_calls/{sample}/{chrom}.vcf"
    log:
        "log/{sample}/call_SNVs_bcftools_chrom.{chrom}.txt"
    params:
        samtools = config["samtools"],
        bcftools = config["bcftools"]
    shell:
        """
        {params.samtools} mpileup -r {wildcards.chrom} -g -f {input.fa} {input.bam} \
        | {params.bcftools} call -mv - | {params.bcftools} view --genotype het --types snps - > {output} 2> {log}
        """

rule merge_SNV_calls:
    input:
        expand("snv_calls/{{sample}}/{chrom}.vcf", chrom = config['chromosomes'])
    output:
        "snv_calls/{sample}/all.vcf"
    shell:
        config["bcftools"] + " concat -O v -o {output} {input}"

rule split_external_snv_calls:
    input:
        vcf = lambda wc: config["snv_calls"][wc.sample],
        tbi = lambda wc: config["snv_calls"][wc.sample] + ".tbi"
    output:
        vcf = "external_snv_calls/{sample}/{chrom}.vcf"
    log: "log/{sample}/external_snv_calls.{chrom}.vcf.log"
    params:
        bcftools = config["bcftools"]
    shell:
        "({params.bcftools} view --samples {wildcards.sample} --types snps {input.vcf} {wildcards.chrom} | {params.bcftools} view --genotype het - > {output.vcf}) > {log} 2>&1"

