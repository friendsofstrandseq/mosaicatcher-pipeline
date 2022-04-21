
################################################################################
# MosaiClassifier                                                              #
################################################################################

rule mosaiClassifier_calc_probs:
    input:
        counts = config["output_location"] + "counts/{sample}/{sample}.txt.gz",
        info   = config["output_location"] + "counts/{sample}/{sample}.info",
        
        states = config["output_location"] + "strandphaser/{sample}/StrandPhaseR_final_output.txt",
        bp     = config["output_location"] + "segmentation/{sample}/Selection_jointseg.txt"
    output:
        output = config["output_location"] + "mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata"
    log:
        config["output_location"] + "log/mosaiClassifier_calc_probs/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier.snakemake.R"

rule create_haplotag_likelihoods:
    input:
        haplotag_table = config["output_location"] + "haplotag/table/{sample}/haplotag_counts_merged.tsv",
        sv_probs_table = config["output_location"] + "mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata",
    output: 
        config["output_location"] + 'mosaiclassifier/haplotag_likelihoods/{sample}.Rdata'
    log:
        config["output_location"] + "log/create_haplotag_likelihoods/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/mosaiclassifier_scripts/haplotagProbs.snakemake.R"

rule mosaiClassifier_make_call:
    input:
        probs = config["output_location"] + 'mosaiclassifier/haplotag_likelihoods/{sample}.Rdata'
    output:
        config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.txt"
    log:
        config["output_location"] + "log/mosaiClassifier_make_call/{sample}/llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.log"
    conda:
        "../envs/rtools.yaml"
    params:
        minFrac_used_bins = 0.8,
        window = config["window"]
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call.snakemake.R"

rule convert_to_html:
    input:
        config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.txt"
    output:
        config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.pdf",
        # report(
        #     config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.html",
        #     category="MosaiClassifier"
        # )
    run:
        pd.read_csv(input[0], sep='\t').to_pdf(output[0])


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
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call_biallelic.snakemake.R"

