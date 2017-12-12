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
        "log/StrandPhaseR.config",
        expand("plots/" + config["sample"] + ".{window}_fixed.pdf", window = [50000, 100000, 200000, 500000]),
        expand("plots/" + config["sample"] + ".{window}_variable.pdf", window = [50000, 100000]),
        expand("segmentation2/" + config["sample"] + ".{window}_fixed.{bpdens}.txt",
               window = [50000, 100000, 200000, 500000], bpdens = ["few","medium","many"]),
        expand("segmentation2/" + config["sample"] + ".{window}_variable.{bpdens}.txt",
               window = [50000, 100000], bpdens = ["few","medium","many"]),
        "strand_states/" + config["sample"] + ".final.txt",
        expand("sv_calls/" + config["sample"] + ".{window}_fixed.{bpdens}.SV_probs.pdf",
               window = [50000, 100000, 200000, 500000], bpdens = ["few","medium","many"]),
        expand("sv_calls/" + config["sample"] + ".{window}_variable.{bpdens}.SV_probs.pdf",
               window = [50000, 100000], bpdens = ["few","medium","many"])



################################################################################
# Plots                                                                        #
################################################################################

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

rule plot_SV_calls:
    input:
        counts = "counts/" + config["sample"] + ".{windows}.txt.gz",
        probs  = "sv_probabilities/" + config["sample"] + ".{windows}.{bpdens}/probabilities.txt"
    output:
        "sv_calls/" + config["sample"] + ".{windows}.{bpdens}.SV_probs.pdf"
    params:
        plot_command = "Rscript " + config["sv_plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.probs} {output}
        touch {output}
        """


################################################################################
# Read counting                                                                #
################################################################################

rule mosaic_count_fixed:
    input:
        bam = expand("bam/{bam}.bam", bam = BAM),
        bai = expand("bam/{bam}.bam.bai", bam = BAM)
    output:
        counts = "counts/" + config["sample"] + ".{window}_fixed.txt.gz",
        info   = "counts/" + config["sample"] + ".{window}_fixed.info"
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
        counts = "counts/" + config["sample"] + ".{window}_variable.txt.gz",
        info   = "counts/" + config["sample"] + ".{window}_variable.info"
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


################################################################################
# SV classification                                                            #
################################################################################

rule install_MaRyam:
    output:
        "utils/R-packages2/MaRyam/R/MaRyam"
    log:
        "log/maryam-install.log"
    shell:
        """
        Rscript  utils/install_maryam.R > {log} 2>&1
        """

rule run_sv_classification:
    input:
        maryam = "utils/R-packages2/MaRyam/R/MaRyam",
        counts = "counts/" + config["sample"] + ".{windows}.txt.gz",
        info   = "counts/" + config["sample"] + ".{windows}.info",
        states = "strand_states/" + config["sample"] + ".final.txt",
        bp     = "segmentation2/" + config["sample"] + ".{windows}.{bpdens}.txt"
    output:
        outdir = "sv_probabilities/" + config["sample"] + ".{windows}.{bpdens}/",
        out1   = "sv_probabilities/" + config["sample"] + ".{windows}.{bpdens}/allSegCellProbs.table"
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
                utils/R-packages2/
        """

rule convert_SVprob_output:
    input:
        probs = "sv_probabilities/" + config["sample"] + ".{windows}.{bpdens}/allSegCellProbs.table",
        info  = "counts/" + config["sample"] + ".{windows}.info"
    output:
        "sv_probabilities/" + config["sample"] + ".{windows}.{bpdens}/probabilities.txt"
    script:
        "utils/helper.convert_svprob_output.R"


################################################################################
# Strand states & phasing                                                      #
################################################################################

rule determine_initial_strand_states:
    input:
        "counts/" + config["sample"] + ".500000_fixed.txt.gz"
    output:
        "strand_states/" + config["sample"] + ".txt"
    params:
        sce_command = "Rscript " + config["sce_script"]
    shell:
        """
        echo "[Note] This is a dirty hack to remove chrX and chrY."
        tmp=$(mktemp)
        {params.sce_command} {input} $tmp
        grep -vP '^Y|^chrY' $tmp | grep -vP '^X|^chrX' > {output}
        """

# Strandphaser needs a different input format which contains the path names to
# the bam files. This rule extracts this information and prepares an input file.
rule convert_strandphaser_input:
    input:
        states = "strand_states/" + config["sample"] + ".txt",
        info   = "counts/" + config["sample"] + ".500000_fixed.info"
    output:
        "strand_states/" + config["sample"] + ".strandphaser_input.txt"
    script:
        "utils/helper.convert_strandphaser_input.R"

rule install_StrandPhaseR:
    output:
        "utils/R-packages/StrandPhaseR/R/StrandPhaseR"
    log:
        "log/strandphaser-install.log"
    shell:
        """
        Rscript  utils/install_strandphaser.R > {log} 2>&1
        """

rule prepare_strandphaser_config:
    input:
        "strand_states/" + config["sample"] + ".txt"
    output:
        "log/StrandPhaseR.config"
    run:
        # Get chromosomes
        chroms = set()
        with open(input[0]) as f:
            skipfirst = True
            for line in f:
                if skipfirst:
                    skipfirst = False
                    continue
                chroms.add(line.split()[0])
        with open(output[0], "w") as f:
            print("[General]",                    file = f)
            print("numCPU           = 1",         file = f)
            print("chromosomes      = c(", ",".join(["'" + chrom + "'" for chrom in chroms]), ")", file = f)
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
            print("exportVCF        = '", config["sample"], ".txt'", sep = "", file = f)
            print("bsGenome         = 'BSgenome.Hsapiens.UCSC.hg19'", file = f)




rule run_strandphaser:
    input:
        mergedbam    = "snv_calls/merged.bam",
        wcregions    = "strand_states/" + config["sample"] + ".strandphaser_input.txt",
        snppositions = "snv_calls/" + config["sample"] + ".vcf",
        configfile   = "log/StrandPhaseR.config",
        strandphaser = "utils/R-packages/StrandPhaseR/R/StrandPhaseR",
        bamfolder    = "bam"
    output:
        "strand_states/" + config["sample"] + ".strandphaser_output.txt"
    log:
        "log/phased_haps.txt.log"
    shell:
        """
        Rscript utils/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                log/StrandPhaseR_analysis \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \
                > {log} 2>&1
        cp log/StrandPhaseR_analysis/Phased/phased_haps.txt {output}
        """


rule convert_strandphaser_output:
    input:
        phased_states  = "strand_states/" + config["sample"] + ".strandphaser_output.txt",
        initial_states = "strand_states/" + config["sample"] + ".txt",
        info           = "counts/" + config["sample"] + ".500000_fixed.info"
    output:
        "strand_states/" + config["sample"] + ".final.txt"
    script:
        "utils/helper.convert_strandphaser_output.R"


################################################################################
# Call SNVs                                                                    #
################################################################################

rule mergeBams:
    input:
        expand("bam/{bam}.bam", bam=BAM)
    output:
        "snv_calls/merged.bam"
    shell:
        config["samtools"] + " merge {output} {input}"

rule indexMergedBam:
    input:
        "snv_calls/merged.bam"
    output:
        "snv_calls/merged.bam.bai"
    shell:
        config["samtools"] + " index {input}"


rule call_SNVs_bcftools_chrom:
    input:
        fa    = config["reference"],
        chrom = "chroms/{chrom}",
        bam   = "snv_calls/merged.bam",
        bai   = "snv_calls/merged.bam.bai"
    output:
        "snv_calls/" + config["sample"] + ".{chrom}.vcf"
    params:
        samtools = config["samtools"],
        bcftools = config["bcftools"]
    shell:
        """
        {params.samtools} mpileup -r {wildcards.chrom} -g -f {input.fa} {input.bam} \
        | {params.bcftools} call -mv - > {output}
        """

# Write one file per chromosome that should be analysed.
rule prepare_chromosomes:
    input:
        "strand_states/" + config["sample"] + ".txt"
    output:
        dynamic("chroms/{chrom}")
    shell:
        """
        tail -n+2 {input} | cut -f1 | sort | uniq | awk '{{print "chroms/" $1}}' | xargs touch
        """

rule merge_SNV_calls:
    input:
        dynamic("snv_calls/" + config["sample"] + ".{chrom}.vcf")
    output:
        expand("snv_calls/" + config["sample"] + ".vcf")
    shell:
        config["bcftools"] + " concat -O v -o {output} {input}"
