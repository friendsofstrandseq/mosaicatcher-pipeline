import math
from collections import defaultdict

configfile: "Snake.config.json"

SAMPLE,BAM = glob_wildcards("bam/{sample}/selected/{bam}.bam")
SAMPLES = sorted(set(SAMPLE))

CELL_PER_SAMPLE= defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)
for sample,bam in zip(SAMPLE,BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace('.sort.mdup',''))

ALLBAMS_PER_SAMPLE = defaultdict(list)
for sample in SAMPLES:
    ALLBAMS_PER_SAMPLE[sample] = glob_wildcards("bam/{}/all/{{bam}}.bam".format(sample)).bam

print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print("  {}:\t{} cells\t {} selected cells".format(s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])))



import os.path

# Current state of the pipeline:
# ==============================
# * count reads in the BAM files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher


METHODS = [
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.005_regfactor6_filterFALSE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.005_regfactor6_filterTRUE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0.005_regfactor6_filterFALSE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0.005_regfactor6_filterTRUE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterFALSE",
    "simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE",
]

BPDENS = [
    "selected_j{}_s{}_scedist{}".format(joint, single, scedist) for joint in [0.1] for single in [0.5] for scedist in [20]
]

singularity: "docker://smei/mosaicatcher-pipeline:v0.1"

localrules:
    all,
    simul,
    simulate_genome,
    add_vafs_to_simulated_genome,
    link_to_simulated_counts,
    link_to_simulated_strand_states,
    generate_exclude_file_1,
    generate_exclude_file_2,
    link_normalized_info_file,
    prepare_segments,
    split_external_snv_calls,
    prepare_strandphaser_config_per_chrom

rule all:
    input:
        expand("plots/{sample}/{window}_fixed.pdf",      sample = SAMPLES, window = [50000, 100000, 200000, 500000]),
        expand("plots/{sample}/{window}_fixed_norm.pdf", sample = SAMPLES, window = [50000, 100000, 200000]),
        expand("sv_calls/{sample}/{window}_fixed_norm.{bpdens}/plots/sv_calls/{method}.{chrom}.pdf",
               sample = SAMPLE,
               chrom = config["chromosomes"],
               window = [100000],
               bpdens = BPDENS,
               method = METHODS),
        expand("ploidy/{sample}/ploidy.{chrom}.txt", sample = SAMPLES, chrom = config["chromosomes"]),
        expand("sv_calls/{sample}/{window}_fixed_norm.{bpdens}/plots/sv_consistency/{method}.consistency-barplot-{plottype}.pdf",
               sample = SAMPLES,
               window = [100000],
               bpdens = BPDENS,
               method = METHODS,
               plottype = ["byaf","bypos"]),
        expand("sv_calls/{sample}/{window}_fixed_norm.{bpdens}/plots/sv_clustering/{method}-{plottype}.pdf",
               sample = SAMPLES,
               window = [100000],
               bpdens = BPDENS,
               method = METHODS,
               plottype = ["position","chromosome"]),
        expand("halo/{sample}/{window}_{suffix}.json.gz",
               sample = SAMPLES,
               window = [100000],
               suffix = ["fixed", "fixed_norm"]),
        expand("stats-merged/{sample}/stats.tsv", sample = SAMPLES),
        expand("postprocessing/merge/{sample}/{window}_fixed_norm.{bpdens}/{method}.txt",
               sample = SAMPLES,
               window = [100000],
               bpdens = BPDENS,
               method = list(set(m.replace('_filterTRUE','').replace('_filterFALSE','') for m in METHODS))),
#        expand("cell-mixing-eval/{window}_fixed_norm.{bpdens}/{method}.tsv",
#               window = [100000],
#               bpdens = BPDENS,
#               method = METHODS),


################################################################################
# Simulation of count data                                                     #
################################################################################

rule simul:
    input:
        expand("sv_calls/simulation{seed}-{window}/{window}_fixed.{segments}/{method}.{chrom}.pdf",
                seed   = list(range(7)),
                window = [50000],
                segments = ["few","medium"],
                method = METHODS,
                chrom = config["chromosomes"]),
        expand("plots/simulation{seed}-{window}/{window}_fixed.pdf",
                seed   = list(range(7)),
                window = [50000])

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
        "log/simulate_counts/genome{seed}-{window_size}.log"
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
# Ploidy estimation                                                            #
################################################################################

rule estimate_ploidy:
    input:
        "counts/{sample}/100000_fixed.txt.gz"
    output:
        "ploidy/{sample}/ploidy.{chrom}.txt"
    log:
        "log/estimate_ploidy/{sample}/{chrom}.log"
    shell:
        """
        python utils/ploidy-estimator.py --chromosome {wildcards.chrom} {input} > {output} 2> {log}
        """



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
        "log/plot_mosaic_counts/{sample}/{file_name}.log"
    params:
        plot_command = "Rscript " + config["plot_script"]
    shell:
        """
        {params.plot_command} {input.counts} {input.info} {output} > {log} 2>&1
        """

ruleorder: plot_SV_calls_simulated > plot_SV_calls

rule plot_SV_calls:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        calls  = "sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.txt",
        complex = "sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.complex.tsv",
        strand = "strand_states/{sample}/{windows}.{bpdens}/final.txt",
        segments = "segmentation2/{sample}/{windows}.{bpdens}.txt",
        scsegments = "segmentation-singlecell/{sample}/{windows}.{bpdens}.txt",
        grouptrack = "postprocessing/group-table/{sample}/{windows}.{bpdens}/{method}.tsv",
    output:
        "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_calls/{method}_filter{filter,(TRUE|FALSE)}.{chrom}.pdf"
    log:
        "log/plot_SV_calls/{sample}/{windows}.{bpdens}.{method}_filter{filter}.{chrom}.log"
    shell:
        """
        Rscript utils/plot-sv-calls.R \
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

rule plot_SV_calls_simulated:
    input:
        counts = "counts/simulation{seed}-{window}/{window}_fixed.txt.gz",
        calls  = "sv_calls/simulation{seed}-{window}/{window}_fixed.{bpdens}/{method}.txt",
        strand = "strand_states/simulation{seed}-{window}/final.txt",
        segments = "segmentation2/simulation{seed}-{window}/{window}_fixed.{bpdens}.txt",
        truth  = "simulation/variants/genome{seed}-{window}.txt"
    output:
        "sv_calls/simulation{seed}-{window}/{window}_fixed.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_calls/{method}.{chrom}.pdf"
    log:
        "log/plot_SV_calls_simulated/simulation{seed}-{window}/{window}_fixed.{bpdens}.{method}.{chrom}.log"
    shell:
        """
        Rscript utils/plot-sv-calls.R \
            segments={input.segments} \
            strand={input.strand} \
            truth={input.truth} \
            calls={input.calls} \
            {input.counts} \
            {wildcards.chrom} \
            {output} 2>&1 > {log}
        """


rule plot_SV_consistency_barplot:
    input:
        sv_calls  = "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt",
    output:
        barplot_bypos = "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_consistency/{method}.consistency-barplot-bypos.pdf",
        barplot_byaf = "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_consistency/{method}.consistency-barplot-byaf.pdf",
    log:
        "log/plot_SV_consistency/{sample}/{windows}.{bpdens}.{method}.log"
    script:
        "utils/sv_consistency_barplot.snakemake.R"


rule generate_halo_json:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
    output:
        json = "halo/{sample}/{windows}.json.gz",
    log:
        "log/generate_halo_json/{sample}/{windows}.{windows}.log"
    shell:
        "(./utils/counts_to_json.py {input.counts} | gzip > {output.json}) 2> {log}"



rule plot_clustering:
    input:
        sv_calls  = "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt",
        binbed = "utils/bin_200kb_all.bed",
    output:
        position = "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_clustering/{method}-position.pdf",
        chromosome = "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/plots/sv_clustering/{method}-chromosome.pdf",
    log:
        "log/plot_clustering/{sample}/{windows}.{bpdens}.{method}.log"
    script:
        "utils/plot-clustering.snakemake.R"


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
        {params.samtools} view -H {input.bam} | awk '/^@SQ/' > {output} 2> {log}
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
                    if contig not in params.chroms:
                        print(contig, file = out)


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
        {params.mc_command} count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {wildcards.window} \
            {input.bam} \
        > {log} 2>&1
        """

rule mosaic_count_variable:
    input:
        bam = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam", bam = BAM_PER_SAMPLE[wc.sample]),
        bai = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam.bai", bam = BAM_PER_SAMPLE[wc.sample]),
        bed = lambda wc: config["variable_bins"][str(wc.window)],
        excl = "log/exclude_file"
    output:
        counts = "counts/{sample}/{window}_variable.txt.gz",
        info   = "counts/{sample}/{window}_variable.info"
    log:
        "log/{sample}/mosaic_count_variable.{window}.log"
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

rule extract_single_cell_counts:
    input:
        "counts/{sample}/{window}_{file_name}.txt.gz"
    output:
        "counts-per-cell/{sample}/{cell}/{window,[0-9]+}_{file_name}.txt.gz"
    shell:
        "zcat {input} | awk '(NR==1) || $5 ==\"{wildcards.cell}\"' | gzip > {output}"


################################################################################
# Normalize counts                                                             #
################################################################################

rule merge_blacklist_bins:
    input:
        norm = "utils/normalization/HGSVC.{window}.txt",
        whitelist = "utils/normalization/inversion-whitelist.tsv",
    output:
        merged = "normalizations/HGSVC.{window}.merged.tsv"
    log:
        "log/merge_blacklist_bins/{window}.log"
    shell:
        """
        utils/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2> {log}
        """

rule normalize_counts:
    input:
        counts = "counts/{sample}/{window}_fixed.txt.gz",
        norm   = "normalizations/HGSVC.{window}.merged.tsv",
    output:
        "counts/{sample}/{window}_fixed_norm.txt.gz"
    log:
        "log/normalize_counts/{sample}/{window}_fixed.log"
    shell:
        """
        Rscript utils/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
        """

rule link_normalized_info_file:
    input:
        info = "counts/{sample}/{window}_fixed.info"
    output:
        info = "counts/{sample}/{window}_fixed_norm.info"
    run:
        d = os.path.dirname(output.info)
        file = os.path.basename(output.info)
        shell("cd {d} && ln -s ../../{input.info} {file} && cd ../..")


################################################################################
# Segmentation                                                                 #
################################################################################

rule segmentation:
    input:
        "counts/{sample}/{window}_{file_name}.txt.gz"
    output:
        "segmentation/{sample}/{window,\d+}_{file_name}.txt.fixme"
    log:
        "log/segmentation/{sample}/{window}_{file_name}.log"
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

# TODO: This is a workaround because latest versions of "mosaic segment" don't compute the "bps"
# TODO: column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    input:
        "segmentation/{sample}/{window}_{file_name}.txt.fixme"
    output:
        "segmentation/{sample}/{window,\d+}_{file_name}.txt"
    shell:
        'awk \'BEGIN {{OFS="\\t"}} {{if ($1=="{wildcards.sample}") $12=int(($14-1)/100000); print}}\' {input} > {output}'

# Pick a few segmentations and prepare the input files for SV classification
rule prepare_segments:
    input:
        "segmentation/{sample}/{windows}.txt"
    output:
        "segmentation2/{sample}/{windows}.{bpdens,(many|medium|few)}.txt"
    log:
        "log/prepare_segments/{sample}/{windows}.{bpdens}.log"
    params:
        quantile = lambda wc: config["bp_density"][wc.bpdens]
    script:
        "utils/helper.prepare_segments.R"

rule segment_one_cell:
    input:
        "counts-per-cell/{sample}/{cell}/{window}_{file_name}.txt.gz"
    output:
        "segmentation-per-cell/{sample}/{cell}/{window,\d+}_{file_name}.txt.fixme"
    log:
        "log/segmentation-per-cell/{sample}/{cell}/{window}_{file_name}.log"
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

# TODO: This is a workaround because latest versions of "mosaic segment" don't compute the "bps"
# TODO: column properly. Remove once fixed in the C++ code.
rule fix_segmentation_one_cell:
    input:
        "segmentation-per-cell/{sample}/{cell}/{window}_{file_name}.txt.fixme"
    output:
        "segmentation-per-cell/{sample}/{cell}/{window,\d+}_{file_name}.txt"
    shell:
        'awk \'BEGIN {{OFS="\\t"}} {{if ($1=="{wildcards.sample}") $12=int(($14-1)/100000); print}}\' {input} > {output}'

rule segmentation_selection:
    input:
        counts="counts/{sample}/{window}_{file_name}.txt.gz",
        jointseg="segmentation/{sample}/{window}_{file_name}.txt",
        singleseg=lambda wc: ["segmentation-per-cell/{}/{}/{}_{}.txt".format(wc.sample, cell, wc.window, wc.file_name) for cell in CELL_PER_SAMPLE[wc.sample]],
        info="counts/{sample}/{window}_{file_name}.info",
    output:
        jointseg="segmentation2/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        singleseg="segmentation-singlecell/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.txt",
        strand_states="strand_states/{sample}/{window,[0-9]+}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}/intitial_strand_state",
    log:
        "log/segmentation_selection/{sample}/{window}_{file_name}.selected_j{min_diff_jointseg}_s{min_diff_singleseg}_scedist{additional_sce_cutoff}.log"
    params:
        cellnames = lambda wc: ",".join(cell for cell in CELL_PER_SAMPLE[wc.sample]),
        sce_min_distance = 500000,
    shell:
        "./utils/detect_strand_states.py --sce_min_distance {params.sce_min_distance} --sce_add_cutoff {wildcards.additional_sce_cutoff}000000 --min_diff_jointseg {wildcards.min_diff_jointseg} --min_diff_singleseg {wildcards.min_diff_singleseg} --output_jointseg {output.jointseg} --output_singleseg {output.singleseg} --output_strand_states {output.strand_states} --samplename {wildcards.sample} --cellnames {params.cellnames} {input.info} {input.counts} {input.jointseg} {input.singleseg} > {log} 2>&1"


################################################################################
# SV classification                                                            #
################################################################################

rule plot_heatmap:
    input:
        maryam = "utils/R-packages2/MaRyam/R/MaRyam",
        haplotypeProbs = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellProbs.table",
        genotypeProbs  = "sv_probabilities/{sample}/{windows}.{bpdens}/allSegCellGTprobs.table",
        info     = "counts/{sample}/{windows}.info",
        bamNames = "sv_probabilities/{sample}/{windows}.{bpdens}/bamNames.txt"
    output:
        "sv_probabilities/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/final_plots/heatmapPlots.pdf"
    params:
        r_package_path = "utils/R-packages2"
    log:
        "log/plot_heatmap/{sample}/{windows}.{bpdens}.log"
    script:
        "utils/plot_heatmap.R"


################################################################################
# New SV classification based on a combination of Sascha's and Maryam's method #
################################################################################

rule mosaiClassifier_make_call:
    input:
        probs = 'haplotag/table/{sample}/haplotag-likelihoods.{window}_fixed_norm.{bpdens}.Rdata'
    output:
        "sv_calls/{sample}/{window}_fixed_norm.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filterFALSE.txt"
    params:
        minFrac_used_bins = 0.8
    log:
        "log/mosaiClassifier_make_call/{sample}/{window}_fixed_norm.{bpdens}.llr{llr}.poppriors{pop_priors}.haplotags{use_haplotags}.gtcutoff{gtcutoff}.regfactor{regfactor}.log"
    script:
        "utils/mosaiClassifier_call.snakemake.R"

rule filter_calls:
    input: 
        calls = "sv_calls/{sample}/{window}_fixed_norm.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.txt"
    output: 
        calls = "sv_calls/{sample}/{window}_fixed_norm.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filterTRUE.txt"
    shell:
        'utils/filter_MosaiCatcher_calls.pl {input.calls} | awk \'BEGIN {{OFS="\\t"}} (NR==1) || ($16=="PASS") {{$16=""; print}}\' > {output.calls}'


rule mosaiClassifier_calc_probs:
    input:
        counts = "counts/{sample}/{windows}.txt.gz",
        info   = "counts/{sample}/{windows}.info",
        states = "strand_states/{sample}/{windows}.{bpdens}/final.txt",
        bp     = "segmentation2/{sample}/{windows}.{bpdens}.txt"
    output:
        output = "sv_probabilities/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/probabilities.Rdata"
    log:
        "log/mosaiClassifier_calc_probs/{sample}/{windows}.{bpdens}.log"
    script:
        "utils/mosaiClassifier.snakemake.R"

rule mosaiClassifier_make_call_biallelic:
    input:
        probs = "sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata"
    output:
        "sv_calls/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/biAllelic_llr{llr}.txt"
    log:
        "log/mosaiClassifier_make_call_biallelic/{sample}/{windows}.{bpdens}.{llr}.log"
    script:
        "utils/mosaiClassifier_call_biallelic.snakemake.R"

rule call_complex_regions:
    input:
        calls  = "sv_calls/{sample}/{windows}.{bpdens}/{method}.txt",
    output:
        complex = "sv_calls/{sample}/{windows}.{bpdens}/{method}.complex.tsv",
    log:
        "log/call_complex_regions/{sample}/{windows}.{bpdens}.{method}.log"
    shell:
        "utils/call-complex-regions.py --merge_distance 5000000 --ignore_haplotypes --min_cell_count 2 {input.calls} > {output.complex} 2>{log}"


rule postprocessing_filter:
    input: 
        calls = "sv_calls/{sample}/{window}_fixed_norm.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.txt"
    output: 
        calls = "postprocessing/filter/{sample}/{window}_fixed_norm.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.txt"
    shell:
        'utils/filter_MosaiCatcher_calls.pl {input.calls}  > {output.calls}'

rule postprocessing_merge:
    input: 
        calls = "postprocessing/filter/{sample}/{window}_fixed_norm.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.txt"
    output: 
        calls = "postprocessing/merge/{sample}/{window}_fixed_norm.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.txt"
    shell:
        'utils/group_nearby_calls_of_same_AF_and_generate_output_table.pl {input.calls}  > {output.calls}'


rule postprocessing_sv_group_table:
    input: 
        calls = "postprocessing/merge/{sample}/{window}_fixed_norm.{bpdens}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.txt"
    output: 
        grouptrack = "postprocessing/group-table/{sample}/{window}_fixed_norm.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}.tsv"
    shell:
        'utils/create-sv-group-track.py {input.calls}  > {output.grouptrack}'


################################################################################
# Strand states & phasing                                                      #
################################################################################

# DEPRECATED rule for calling SCEs / determining the initial strand states using 
# the C++ MosaiCatcher code.
#rule determine_initial_strand_states:
    #input:
        #"counts/{sample}/500000_fixed.txt.gz"
    #output:
        #"strand_states/{sample}/intitial_strand_state"
    #log:
        #"log/determine_initial_strand_states/{sample}.log"
    #params:
        #mc_command = config["mosaicatcher"]
    #shell:
        #"""
        #{params.mc_command} states -o {output} {input} 2>&1 > {log}
        #"""

#rule determine_initial_strand_states:
    #input:
        #"strand_states/{sample}/100000_fixed_norm.selected_j0.5_s1.0.intitial_strand_state"
    #output:
        #"strand_states/{sample}/intitial_strand_state"
    #shell:
        #"""
        #cd strand_states/{wildcards.sample} && ln -s 100000_fixed_norm.selected_j0.5_s1.0.intitial_strand_state intitial_strand_state && cd ../..
        #"""

# Strandphaser needs a different input format which contains the path names to
# the bam files. This rule extracts this information and prepares an input file.
rule convert_strandphaser_input:
    input:
        states = "strand_states/{sample}/{windows}.{bpdens}/intitial_strand_state",
        info   = "counts/{sample}/500000_fixed.info"
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/strandphaser_input.txt"
    log:
        "log/convert_strandphaser_input/{sample}/{windows}.{bpdens}.log"
    script:
        "utils/helper.convert_strandphaser_input.R"

rule install_StrandPhaseR:
    output:
        "utils/R-packages/StrandPhaseR/R/StrandPhaseR"
    log:
        "log/install_StrandPhaseR.log"
    shell:
        """
        TAR=$(which tar) Rscript utils/install_strandphaser.R > {log} 2>&1
        """

rule prepare_strandphaser_config_per_chrom:
    input:
        "strand_states/{sample}/{windows}.{bpdens}/intitial_strand_state"
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR.{chrom}.config"
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


def locate_snv_vcf(wildcards):
    if "snv_calls" not in config or wildcards.sample not in config["snv_calls"] or config["snv_calls"][wildcards.sample] == "":
        if "snv_sites_to_genotype" in config and config["snv_sites_to_genotype"] != "":
            return "snv_genotyping/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
        else:
            return "snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)
    else:
        return "external_snv_calls/{}/{}.vcf".format(wildcards.sample, wildcards.chrom)

rule run_strandphaser_per_chrom:
    input:
        wcregions    = "strand_states/{sample}/{windows}.{bpdens}/strandphaser_input.txt",
        snppositions = locate_snv_vcf,
        configfile   = "strand_states/{sample}/{windows}.{bpdens}/StrandPhaseR.{chrom}.config",
        strandphaser = "utils/R-packages/StrandPhaseR/R/StrandPhaseR",
        bamfolder    = "bam/{sample}/selected"
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf",
    log:
        "log/run_strandphaser_per_chrom/{sample}/{windows}.{bpdens}/{chrom}.log"
    shell:
        """
        Rscript utils/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                strand_states/{wildcards.sample}/{wildcards.windows}.{wildcards.bpdens}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \
                > {log} 2>&1
        """

rule compress_vcf:
    input:
        vcf="{file}.vcf",
    output:
        vcf="{file}.vcf.gz",
    log:
        "log/compress_vcf/{file}.log"
    shell:
        "(cat {input.vcf} | bgzip > {output.vcf}) > {log} 2>&1"


rule index_vcf:
    input:
        vcf="{file}.vcf.gz",
    output:
        tbi="{file}.vcf.gz.tbi",
    shell:
        "bcftools index --tbi {input.vcf}"

rule merge_strandphaser_vcfs:
    input:
        vcfs=expand("strand_states/{{sample}}/{{windows}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"]),
        tbis=expand("strand_states/{{sample}}/{{windows}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"]),
    output:
        vcf='phased-snvs/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.vcf.gz'
    log:
        "log/merge_strandphaser_vcfs/{sample}/{windows}.{bpdens}.log"
    shell:
        "(bcftools concat -a {input.vcfs} | bcftools view -o {output.vcf} -O z --genotype het --types snps - ) > {log} 2>&1"



rule combine_strandphaser_output:
    input:
        expand("strand_states/{{sample}}/{{windows}}.{{bpdens}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
                chrom = config["chromosomes"])
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/strandphaser_output.txt"
    log:
        "log/combine_strandphaser_output/{sample}/{windows}.{bpdens}.log"
    shell:
        """
        set +o pipefail
        cat {input} | head -n1 > {output};
		tail -q -n+2 {input} >> {output};
        """


rule convert_strandphaser_output:
    input:
        phased_states  = "strand_states/{sample}/{windows}.{bpdens}/strandphaser_output.txt",
        initial_states = "strand_states/{sample}/{windows}.{bpdens}/intitial_strand_state",
        info           = "counts/{sample}/500000_fixed.info"
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/final.txt"
    log:
        "log/convert_strandphaser_output/{sample}/{windows}.{bpdens}.log"
    script:
        "utils/helper.convert_strandphaser_output.R"


################################################################################
# Haplotagging                                                                 #
################################################################################

rule haplotag_bams:
    input:
        vcf='phased-snvs/{sample}/{windows}.{bpdens}.vcf.gz',
        tbi='phased-snvs/{sample}/{windows}.{bpdens}.vcf.gz.tbi',
        bam='bam/{sample}/selected/{bam}.bam',
        bai='bam/{sample}/selected/{bam}.bam.bai',
        ref = config["reference"],
    output:
        bam='haplotag/bam/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/{bam}.bam',
    log:
        "log/haplotag_bams/{sample}/{windows}.{bpdens}/{bam}.log"
    shell:
        "whatshap haplotag -o {output.bam} -r {input.ref} {input.vcf} {input.bam} > {log} 2>{log}"

rule create_haplotag_segment_bed:
    input:
        segments="segmentation2/{sample}/{size}{what}.{bpdens}.txt",
    output:
        bed="haplotag/bed/{sample}/{size,[0-9]+}{what}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.bed",
    shell:
        "awk 'BEGIN {{s={wildcards.size};OFS=\"\\t\"}} $2!=c {{prev=0}} NR>1 {{print $2,prev*s+1,($3+1)*s; prev=$3+1; c=$2}}' {input.segments} > {output.bed}"

rule create_haplotag_table:
    input:
        bam='haplotag/bam/{sample}/{windows}.{bpdens}/{cell}.bam',
        bai='haplotag/bam/{sample}/{windows}.{bpdens}/{cell}.bam.bai',
        bed = "haplotag/bed/{sample}/{windows}.{bpdens}.bed"
    output:
        tsv='haplotag/table/{sample}/by-cell/haplotag-counts.{cell}.{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.tsv'
    log:
        "log/create_haplotag_table/{sample}.{cell}.{windows}.{bpdens}.log"
    script:
        "utils/haplotagTable.snakemake.R"

rule merge_haplotag_tables:
    input:
        tsvs=lambda wc: ['haplotag/table/{}/by-cell/haplotag-counts.{}.{}.{}.tsv'.format(wc.sample,cell,wc.windows,wc.bpdens) for cell in BAM_PER_SAMPLE[wc.sample]],
    output:
        tsv='haplotag/table/{sample}/full/haplotag-counts.{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.tsv'
    shell:
        '(head -n1 {input.tsvs[0]} && tail -q -n +2 {input.tsvs}) > {output.tsv}'


rule create_haplotag_likelihoods:
    input:
        haplotag_table='haplotag/table/{sample}/full/haplotag-counts.{windows}.{bpdens}.tsv',
        sv_probs_table = 'sv_probabilities/{sample}/{windows}.{bpdens}/probabilities.Rdata',
    output: 'haplotag/table/{sample}/haplotag-likelihoods.{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.Rdata'
    log:
        "log/create_haplotag_likelihoods/{sample}.{windows}.{bpdens}.log"
    script:
        "utils/haplotagProbs.snakemake.R"


################################################################################
# Call SNVs                                                                    #
################################################################################

rule mergeBams:
    input:
        lambda wc: expand("bam/" + wc.sample + "/all/{bam}.bam", bam = ALLBAMS_PER_SAMPLE[wc.sample]) if wc.sample in ALLBAMS_PER_SAMPLE else "FOOBAR",
    output:
        "snv_calls/{sample}/merged.bam"
    log:
        "log/mergeBams/{sample}.log"
    threads:
        4
    shell:
        config["samtools"] + " merge -@ {threads} {output} {input} 2>&1 > {log}"

rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    log:
        "{file}.bam.log"
    shell:
        config["samtools"] + " index {input} 2> {log}"

rule call_SNVs_bcftools_chrom:
    input:
        bam   = "snv_calls/{sample}/merged.bam",
        bai   = "snv_calls/{sample}/merged.bam.bai"
    output:
        "snv_calls/{sample}/{chrom,chr[0-9A-Z]+}.vcf"
    log:
        "log/call_SNVs_bcftools_chrom/{sample}/{chrom}.log"
    params:
        fa = config["reference"],
        samtools = config["samtools"],
        bcftools = config["bcftools"]
    shell:
        """
        {params.samtools} mpileup -r {wildcards.chrom} -g -f {params.fa} {input.bam} \
        | {params.bcftools} call -mv - | {params.bcftools} view --genotype het --types snps - > {output} 2> {log}
        """

rule regenotype_SNVs:
    input:
        bam   = "snv_calls/{sample}/merged.bam",
        bai   = "snv_calls/{sample}/merged.bam.bai",
        fa = config["reference"],
        sites = config["snv_sites_to_genotype"],
    output:
        vcf = "snv_genotyping/{sample}/{chrom,chr[0-9A-Z]+}.vcf"
    log:
        "log/snv_genotyping/{sample}/{chrom}.log"
    params:
        bcftools = config["bcftools"]
    shell:
        """
        (freebayes -f {input.fa} -r {wildcards.chrom} -@ {input.sites} --only-use-input-alleles {input.bam} --genotype-qualities | {params.bcftools} view --exclude-uncalled --genotype het --types snps --include "QUAL>=10" - > {output.vcf}) 2> {log}
        """

rule merge_SNV_calls:
    input:
        expand("snv_calls/{{sample}}/{chrom}.vcf", chrom = config['chromosomes'])
    output:
        "snv_calls/{sample}/all.vcf"
    log:
        "log/merge_SNV_calls/{sample}.log"
    shell:
        config["bcftools"] + " concat -O v -o {output} {input} 2>&1 > {log}"

rule split_external_snv_calls:
    input:
        vcf = lambda wc: config["snv_calls"][wc.sample],
        tbi = lambda wc: config["snv_calls"][wc.sample] + ".tbi"
    output:
        vcf = "external_snv_calls/{sample}/{chrom}.vcf"
    log:
        "log/split_external_snv_calls/{sample}/{chrom}.vcf.log"
    params:
        bcftools = config["bcftools"]
    shell:
        """
        ({params.bcftools} view --samples {wildcards.sample} \
            --types snps \
            --exclude-uncalled \
            --trim-alt-alleles \
            -m 2 -M 2 \
            {input.vcf} \
            {wildcards.chrom} \
        | {params.bcftools} view --genotype het - \
        > {output.vcf} ) \
        > {log} 2>&1
        """


################################################################################
# Summary statistics on sv calls                                               #
################################################################################


rule summary_statistics:
    input:
        segmentation = 'segmentation2/{sample}/{windows}.{bpdens}.txt',
        strandstates = 'strand_states/{sample}/{windows}.{bpdens}/intitial_strand_state',
        sv_calls = 'sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.txt',
        complex = "sv_calls/{sample}/{windows}.{bpdens}/{method}_filter{filter}.complex.tsv",
        merged = "postprocessing/merge/{sample}/{windows}.{bpdens}/{method}.txt",
    output:
        tsv = 'stats/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/{method}_filter{filter,(TRUE|FALSE)}.tsv',
    log:
        'log/summary_statistics/{sample}/{windows}.{bpdens}/{method}_filter{filter}.log'
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
        shell('utils/callset_summary_stats.py --segmentation {input.segmentation} --strandstates {input.strandstates} --complex-regions {input.complex} {additional_params} {input.sv_calls}  > {output.tsv} 2> {log}')

rule aggregate_summary_statistics:
    input:
        tsv=expand("stats/{{sample}}/{window}_fixed_norm.{bpdens}/{method}.tsv", window = [100000], bpdens = BPDENS, method = METHODS),
    output:
        tsv="stats-merged/{sample}/stats.tsv"
    shell:
        "(head -n1 {input.tsv[0]} && (tail -n1 -q {input.tsv} | sort -k1) ) > {output}"
    
rule evaluate_cell_mixing:
    input:
        sv_calls = expand('sv_calls/{sample}/{{windows}}.{{bpdens}}/{{method}}.txt', sample=SAMPLES),
        truth = '../input-data/ground_truth/RPE-BM510_manual/clonal-events.tsv',
    output:
        tsv = 'cell-mixing-eval/{windows}.{bpdens}/{method}.tsv'
    log:
        'cell-mixing-eval/{windows}.{bpdens}/{method}.log'
    run:
        names = ','.join(SAMPLES)
        shell('utils/evaluate_cell_mixing.py --names {names} {input.truth} {input.sv_calls} > {output.tsv} 2> {log}')
