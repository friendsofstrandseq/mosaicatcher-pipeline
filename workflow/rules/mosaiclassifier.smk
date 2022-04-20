
################################################################################
# MosaiClassifier                                                              #
################################################################################

rule mosaiClassifier_calc_probs:
    input:
        counts = config["output_location"] + "counts/{sample}/{sample}.txt.gz",
        info   = config["output_location"] + "counts/{sample}/{sample}.info",
        
        states = config["output_location"] + "strandphaser/{sample}/StrandPhaseR_final_output.txt",
        bp     = config["output_location"] + "segmentation/{sample}/{sample}.txt"
    output:
        output = config["output_location"] + "mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata"
    log:
        config["output_location"] + "log/mosaiClassifier_calc_probs/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "utils/mosaiClassifier.snakemake.R"

rule create_haplotag_likelihoods:
    input:
        haplotag_table = config["output_location"] + "haplotag/table/{sample}/haplotag_counts_merged.tsv",
        sv_probs_table = config["output_location"] + "mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata",
    output: 
        config["output_location"] + 'mosaiclassifier/table/{sample}/haplotag-likelihoods.{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}.Rdata'
    log:
        config["output_location"] + "log/create_haplotag_likelihoods/{sample}.{window}.{bpdens}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "utils/haplotagProbs.snakemake.R"

rule mosaiClassifier_make_call:
    input:
        probs = config["output_location"] + 'haplotag/table/{sample}/haplotag-likelihoods.{window}.{bpdens}.Rdata'
    output:
        config["output_location"] + "sv_calls/{sample}/{window}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/simpleCalls_llr{llr}_poppriors{pop_priors,(TRUE|FALSE)}_haplotags{use_haplotags,(TRUE|FALSE)}_gtcutoff{gtcutoff,[0-9\\.]+}_regfactor{regfactor,[0-9]+}_filterFALSE.txt"
    log:
        config["output_location"] + "log/mosaiClassifier_make_call/{sample}/{window}.{bpdens}.llr{llr}.poppriors{pop_priors}.haplotags{use_haplotags}.gtcutoff{gtcutoff}.regfactor{regfactor}.log"
    conda:
        "../envs/rtools.yaml"
    params:
        minFrac_used_bins = 0.8
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
    conda:
        "../envs/rtools.yaml"
    script:
        "utils/mosaiClassifier_call_biallelic.snakemake.R"

