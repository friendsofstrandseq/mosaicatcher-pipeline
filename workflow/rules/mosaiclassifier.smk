# from workflow.scripts.utils.utils import get_mem_mb

################################################################################
# MosaiClassifier                                                              #
################################################################################


rule mosaiClassifier_calc_probs:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        info="{folder}/{sample}/counts/{sample}.info",
        states="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
        bp="{folder}/{sample}/segmentation/Selection_jointseg.txt",
    output:
        output="{folder}/{sample}/mosaiclassifier/sv_probabilities/probabilities.Rdata",
    log:
        "{folder}/log/mosaiClassifier_calc_probs/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier.snakemake.R"


rule create_haplotag_likelihoods:
    input:
        haplotag_table="{folder}/{sample}/haplotag/table/haplotag_counts_merged.tsv",
        sv_probs_table="{folder}/{sample}/mosaiclassifier/sv_probabilities/probabilities.Rdata",
    output:
        "{folder}/{sample}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
    log:
        "{folder}/log/create_haplotag_likelihoods/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/haplotagProbs.snakemake.R"


# rule mosaiClassifier_make_call:
#     input:
#         probs="{folder}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
#     output:
#         "{folder}/mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
#     log:
#         "{folder}/log/mosaiClassifier_make_call/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log",
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
        probs="{folder}/{sample}/mosaiclassifier/haplotag_likelihoods/{sample}.Rdata",
    output:
        "{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filterFALSE.tsv",
    log:
        "{folder}/log/mosaiClassifier_make_call/{sample}/{method}.log",
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
        probs="{folder}/{sample}/sv_probabilities/*probabilities.Rdata",
    output:
        "{folder}/{sample}/sv_calls/biAllelic_llr{llr}.txt",
    log:
        "{folder}/log/mosaiClassifier_make_call_biallelic/{sample}/{llr}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/mosaiclassifier_scripts/mosaiClassifier_call_biallelic.snakemake.R"


# rule debug_complex:
#     input:
#         calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
#     output:
#         "{folder}/{sample}/mosaiclassifier/sv_calls_debug/{method}_filter{filter}.ok"
#     log:
#         "{folder}/log/mosaiclassifier/sv_calls_debug/{sample}/{method}_filter{filter}.log",
#     shell:
#         "cat {input.calls}"


rule call_complex_regions:
    input:
        calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        # debug="{folder}/{sample}/mosaiclassifier/sv_calls_debug/{method}_filter{filter}.ok",
    output:
        complex_regions="{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv",
    log:
        "{folder}/log/call_complex_regions/{sample}/{method}_filter{filter}.log",
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
