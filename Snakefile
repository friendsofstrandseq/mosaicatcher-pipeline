configfile: "Snake.config.json"

SAMPLE,BAM = glob_wildcards("bam/{sample}/{bam}.bam")
BAM_PER_SAMPLE = dict([(s,[]) for s in SAMPLE])
for i in range(len(SAMPLE)):
    BAM_PER_SAMPLE[SAMPLE[i]].append(BAM[i])
print("Detected {} samples:".format(len(set(SAMPLE))))
for s in set(SAMPLE):
    print("  {}:\t{} cells".format(s, len(BAM_PER_SAMPLE[s])))

import os.path

# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher


METHODS = ["simpleCalls"]

rule all:
    input:
        expand("plots/{sample}/{window}_fixed.pdf", sample = SAMPLE, window = [50000, 100000, 200000, 500000]),
        expand("plots/{sample}/{window}_variable.pdf", sample = SAMPLE, window = [50000, 100000]),
        expand("sv_calls/{sample}/{window}_fixed.{bpdens}/{method}.{chrom}.pdf",
               sample = SAMPLE,
               chrom = config["chromosomes"],
               window = [100000],
               bpdens = ["few","medium","many"],
               method = METHODS)


################################################################################
# Simulation of count data                                                     #
################################################################################

rule simul:
    input:
        expand("sv_calls/simulation{seed}-{window}/{window}_fixed.{segments}/{method}.{chrom}.pdf",
                seed   = list(range(7)),
                window = [50000],
                segments = ["medium"],
                method = METHODS,
                chrom = config["chromosomes"]),
        expand("plots/simulation{seed}-{window}/{window}_fixed.pdf",
                seed   = list(range(7)),
                window = [50000]),
        expand("evaluation/simulation{seed}_{window}.{segments}.pdf",
                seed  = list(range(7)),
                window = [50000],
                segments = ["medium"]),
        expand("sv_probabilities/simulation{seed}-{window}/{window}_fixed.{segments}/final_plots/heatmapPlots.pdf",
                seed  = list(range(7)),
                window = [50000],
                segments = ["medium"])

rule evaluate_simulation:
    input:
        prob = "sv_probabilities/simulation{seed}-{window}/{window}_fixed.{bpdens}/probabilities.txt",
        simul = "simulation/variants/genome{seed}-{window}.txt"
    output:
        "evaluation/simulation{seed}_{window}.{bpdens}.pdf"
    script:
        "utils/evaluate_simulation.R"

rule simulate_genome:
    output:
        tsv="simulation/genome/genome{seed}.tsv"
    log:
        "log/simulate_genome/genome{seed}.tsv"
    params:
        svcount     =     200,
        minsize     =  100000,
        maxsize     = 5000000,
        mindistance = 1000000,
    shell:
        "utils/simulate_SVs.R {wildcards.seed} {params.svcount} {params.minsize} {params.maxsize} {params.mindistance} {output.tsv} > {log} 2>&1"

rule add_vafs_to_simulated_genome:
    input:
        tsv="simulation/genome/genome{seed}.tsv"
    output:
        tsv="simulation/genome-with-vafs/genome{seed}.tsv"
    params:
        min_vaf = config["simulation_min_vaf"],
        max_vaf = config["simulation_max_vaf"],
    shell:
        """
        awk -v min_vaf={params.min_vaf} -v max_vaf={params.max_vaf} -v seed={wildcards.seed} \
        'BEGIN {{srand(seed); OFS="\\t"}} {{vaf=min_vaf+rand()*(max_vaf-min_vaf); print $0, vaf}}' {input.tsv} > {output.tsv}
        """

def min_coverage(wildcards):
    return round(float(config["simulation_min_reads_per_library"]) * int(wildcards.window_size) / float(config["genome_size"]))

def max_coverage(wildcards):
    return round(float(config["simulation_max_reads_per_library"]) * int(wildcards.window_size) / float(config["genome_size"]))

def neg_binom_p(wildcards):
    return float(config["simulation_neg_binom_p"][wildcards.window_size])

rule simulate_counts:
    input:
        config="simulation/genome-with-vafs/genome{seed}.tsv",
    output:
        counts="simulation/counts/genome{seed}-{window_size}.txt.gz",
        segments="simulation/segments/genome{seed}-{window_size}.txt",
        phases="simulation/phases/genome{seed}-{window_size}.txt",
        info="simulation/info/genome{seed}-{window_size}.txt",
        sce="simulation/sce/genome{seed}-{window_size}.txt",
        variants="simulation/variants/genome{seed}-{window_size}.txt",
    params:
        mc_command   = config["mosaicatcher"],
        neg_binom_p  = neg_binom_p,
        min_coverage = min_coverage,
        max_coverage = max_coverage,
        cell_count   = config["simulation_cell_count"],
        alpha        = config["simulation_alpha"],
    log:
        "log/simulate_counts/genome{seed}-{window_size}.txt"
    shell:
        """
            {params.mc_command} simulate \
            -w {wildcards.window_size} \
            --seed {wildcards.seed} \
            -n {params.cell_count} \
            -p {params.neg_binom_p} \
            -c {params.min_coverage} \
            -C {params.max_coverage} \
            -a {params.alpha} \
            -V {output.variants} \
            -i {output.info} \
            -o {output.counts} \
            -U {output.segments} \
            -P {output.phases} \
            -S {output.sce} \
            --sample-name simulation{wildcards.seed}-{wildcards.window_size} \
            {input.config} > {log} 2>&1
        """

rule link_to_simulated_counts:
    input:
        counts="simulation/counts/genome{seed}-{window_size}.txt.gz",
        info="simulation/info/genome{seed}-{window_size}.txt",
    output:
        counts = "counts/simulation{seed}-{window_size}/{window_size}_fixed.txt.gz",
        info   = "counts/simulation{seed}-{window_size}/{window_size}_fixed.info"
    run:
        d = os.path.dirname(output.counts)
        count_file = os.path.basename(output.counts)
        info_file = os.path.basename(output.info)
        shell("cd {d} && ln -s ../../{input.counts} {count_file} && ln -s ../../{input.info} {info_file} && cd ../..")


rule link_to_simulated_strand_states:
    input:
        sce="simulation/sce/genome{seed}-{window_size}.txt",
    output:
        states="strand_states/simulation{seed}-{window_size}/final.txt",
    run:
        d = os.path.dirname(output.states)
        f = os.path.basename(output.states)
        shell("cd {d} && ln -s ../../{input.sce} {f} && cd ../..")

ruleorder: link_to_simulated_counts > mosaic_count_fixed
ruleorder: link_to_simulated_strand_states > convert_strandphaser_output

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
        calls  = "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/{method}.{chrom}.pdf"
    shell:
        """
        Rscript utils/chrom.R calls={input.calls} {input.counts} {wildcards.chrom} {output}
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
        bam = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
        bai = lambda wc: expand("bam/" + wc.sample + "/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
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
        TAR=$(which tar) Rscript utils/install_maryam.R > {log} 2>&1
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
        out1   = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table",
        out2   = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellGTprobs.table",
        bamNames = "sv_probabilities/{sample}/{windows}.{bpdens}/bamNames.txt"
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

rule plot_heatmap:
    input:
        maryam = "utils/R-packages2/MaRyam/R/MaRyam",
        haplotypeProbs = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table",
        genotypeProbs  = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellGTprobs.table",
        info     = "counts/{sample}/{windows}.info",
        bamNames = "sv_probabilities/{sample}/{windows}.{bpdens}/bamNames.txt"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens}/final_plots/heatmapPlots.pdf"
    params:
        r_package_path = "utils/R-packages2"
    log:
        "log/{sample}/final.plots.{windows}.{bpdens}.txt"
    script:
        "utils/plot_heatmap.R"

rule convert_SVprob_output:
    input:
        probs    = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table",
        info     = "counts/{sample}/{windows}.info",
        bamNames = "sv_probabilities/{sample}/{windows}.{bpdens}/bamNames.txt"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.txt"
    params:
        sample_name = lambda wc: wc.sample
    log:
        "log/{sample}/convert_SVprob_output.{windows}.{bpdens}.txt"
    script:
        "utils/helper.convert_svprob_output.R"


# rule run_sv_classification_new:
#     input:
#         counts = "counts/{sample}/{windows}.txt.gz",
#         info   = "counts/{sample}/{windows}.info",
#         states = "strand_states/{sample}/final.txt",
#         bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
#     output:
#         "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.txt"
#     script:
#         "utils/sv_classifier.R"


################################################################################
# New SV classification based on a combination of Sascha's and Maryam's method #
################################################################################

rule mosaiClassifier_make_call:
    input:
        "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens}/simpleCalls.txt"
    shell:
        """
        Rscript utils/mosaiClassifier.makeCall.R {input} {output}
        """

rule mosaiClassifier_calc_probs:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        output = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    script:
        "utils/mosaiClassifier.main.R"


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
        TAR=$(which tar) Rscript utils/install_strandphaser.R > {log} 2>&1
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
        "strand_states/{sample}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt"
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
        expand("strand_states/{{sample}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
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
        lambda wc: expand("bam/" + wc.sample + "/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]) if wc.sample in BAM_PER_SAMPLE else "FOOBAR",
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
        "({params.bcftools} view --samples {wildcards.sample} --types snps --exclude-uncalled --trim-alt-alleles -m 2 -M 2 {input.vcf} {wildcards.chrom} | {params.bcftools} view --genotype het - > {output.vcf}) > {log} 2>&1"

