
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
        config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv"
    log:
        config["output_location"] + "log/mosaiClassifier_make_call/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log"
    conda:
        "../envs/rtools.yaml"
    params:
        minFrac_used_bins = 0.8,
        window = config["window"]
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call.snakemake.R"

# rule convert_to_html:
#     input:
#         config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.txt"
#     output:
#         config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.pdf",
#         # report(
#         #     config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filter{filter}.html",
#         #     category="MosaiClassifier"
#         # )
#     run:
#         pd.read_csv(input[0], sep='\t').to_pdf(output[0])


# CHECKME : check if still useful ?
rule mosaiClassifier_make_call_biallelic:
    input:
        probs = config["output_location"] + "sv_probabilities/{sample}/probabilities.Rdata"
    output:
        config["output_location"] + "sv_calls/{sample}/biAllelic_llr{llr}.txt"
    log:
        config["output_location"] + "log/mosaiClassifier_make_call_biallelic/{sample}/{llr}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call_biallelic.snakemake.R"


rule call_complex_regions:
    input:
        calls  = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}.tsv",
    output:
        complex_regions = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}.complex.tsv",
    log:
        config["output_location"] + "log/call_complex_regions/{sample}/{method}.log"
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        scripts/mosaiclassifier_scripts/call-complex-regions.py \
        --merge_distance 5000000 \
        --ignore_haplotypes \
        --min_cell_count 2 {input.calls} > {output.complex_regions} 2>{log}
        """
