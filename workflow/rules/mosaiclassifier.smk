# from workflow.scripts.utils.utils import get_mem_mb

################################################################################
# MosaiClassifier                                                              #
################################################################################


rule mosaiClassifier_calc_probs:
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        info="{output_folder}/counts/{sample}/{sample}.info",
        states="{output_folder}/strandphaser/{sample}/StrandPhaseR_final_output.txt",
        bp="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
    output:
        output="{output_folder}/mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata",
    log:
        "{output_folder}/log/mosaiClassifier_calc_probs/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier.snakemake.R"


rule create_haplotag_likelihoods:
    input:
        haplotag_table="{output_folder}/haplotag/table/{sample}/haplotag_counts_merged.tsv",
        sv_probs_table="{output_folder}/mosaiclassifier/sv_probabilities/{sample}/probabilities.Rdata",
    output:
        "{output_folder}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
    log:
        "{output_folder}/log/create_haplotag_likelihoods/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/haplotagProbs.snakemake.R"


# rule mosaiClassifier_make_call:
#     input:
#         probs="{output_folder}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
#     output:
#         "{output_folder}/mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
#     log:
#         "{output_folder}/log/mosaiClassifier_make_call/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log",
#     conda:
#         "../envs/rtools.yaml"
#     params:
#         minFrac_used_bins=0.8,
#         window=config["window"],
#     resources:
#         mem_mb=get_mem_mb,
#     script:
#         "../scripts/mosaiclassifier_scripts/mosaiClassifier_call.snakemake.R"


rule mosaiClassifier_make_call:
    input:
        probs="{output_folder}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
    output:
        "{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filterFALSE.tsv",
    log:
        "{output_folder}/log/mosaiClassifier_make_call/{sample}/{method}.log",
    conda:
        "../envs/rtools.yaml"
    params:
        minFrac_used_bins=0.8,
        window=config["window"],
        llr=lambda wc: config["methods"][wc.method]["llr"],
        poppriors=lambda wc: config["methods"][wc.method]["poppriors"],
        haplotags=lambda wc: config["methods"][wc.method]["haplotags"],
        gtcutoff=lambda wc: config["methods"][wc.method]["gtcutoff"],
        regfactor=lambda wc: config["methods"][wc.method]["regfactor"],
        filter=lambda wc: config["methods"][wc.method]["filter"],
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call.snakemake.R"


rule mosaiClassifier_make_call_biallelic:
    input:
        probs="{output_folder}/sv_probabilities/{sample}/probabilities.Rdata",
    output:
        "{output_folder}/sv_calls/{sample}/biAllelic_llr{llr}.txt",
    log:
        "{output_folder}/log/mosaiClassifier_make_call_biallelic/{sample}/{llr}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call_biallelic.snakemake.R"

rule debug_complex:
    input:
        calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
    output:
        touch("{output_folder}/mosaiclassifier/sv_calls_debug/{sample}/{method}_filter{filter}.ok")
    log:
        "{output_folder}/mosaiclassifier/sv_calls_debug/{sample}/{method}_filter{filter}.log"
    shell:
        "cat {input.calls}"

rule call_complex_regions:
    input:
        calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
        debug="{output_folder}/mosaiclassifier/sv_calls_debug/{sample}/{method}_filter{filter}.ok"
    output:
        complex_regions="{output_folder}/mosaiclassifier/complex/{sample}/{method}_filter{filter}.tsv",
    log:
        "{output_folder}/log/call_complex_regions/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        workflow/scripts/mosaiclassifier_scripts/call-complex-regions.py \
        --merge_distance 5000000 \
        --ignore_haplotypes \
        --min_cell_count 2 {input.calls} > {output.complex_regions} 2>{log}
        """
