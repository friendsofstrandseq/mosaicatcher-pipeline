# from workflow.scripts.utils.utils import get_mem_mb
################################################################################
# PostProcessing                                                               #
################################################################################


# CHECKME : segdup file only in hg38
rule postprocessing_filter:
    input:
        calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filterFALSE.tsv",
    output:
        calls="{folder}/{sample}/mosaiclassifier/postprocessing/filter/{method}.tsv",
    log:
        "{folder}/log/mosaiclassifier/postprocessing/filter/{sample}/{method}.log",
    conda:
        "../envs/mc_base.yaml"
    params:
        # segdups=config["segdups"],
        segdups=config["references_data"][config["reference"]]["segdups"],
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        export LC_CTYPE=en_US.UTF-8 
        export LC_ALL=en_US.UTF-8 
        workflow/scripts/postprocessing/filter_MosaiCatcher_calls.pl {input.calls} {params.segdups} > {output.calls}
        """


# # CHECKME : segdup file only in hg38
# rule postprocessing_filter:
#     input:
#         calls="{folder}/mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
#     output:
#         calls="{folder}/mosaiclassifier/postprocessing/filter/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#     log:
#         "{folder}/log/mosaiclassifier/postprocessing/filter/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log",
#     conda:
#         "../envs/mc_base.yaml"
#     params:
#         # segdups=config["segdups"],
#         segdups=config["references_data"][config["reference"]]["segdups"],
#     resources:
#         mem_mb=get_mem_mb,
#     shell:
#         """
#         export LC_CTYPE=en_US.UTF-8
#         export LC_ALL=en_US.UTF-8
#         workflow/scripts/postprocessing/filter_MosaiCatcher_calls.pl {input.calls} {params.segdups} > {output.calls}
#         """


rule postprocessing_merge:
    input:
        calls="{folder}/{sample}/mosaiclassifier/postprocessing/filter/{method}.tsv",
    output:
        calls="{folder}/{sample}/mosaiclassifier/postprocessing/merge/{method}.tsv",
    log:
        "{folder}/log/mosaiclassifier/postprocessing/merge/{sample}/{method}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        export LC_CTYPE=en_US.UTF-8 
        export LC_ALL=en_US.UTF-8 
        workflow/scripts/postprocessing/group_nearby_calls_of_same_AF_and_generate_output_table.pl {input.calls}  > {output.calls}
        """


# rule postprocessing_merge:
#     input:
#         calls="{folder}/mosaiclassifier/postprocessing/filter/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#     output:
#         calls="{folder}/mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#     log:
#         "{folder}/log/mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log",
#     conda:
#         "../envs/mc_base.yaml"
#     resources:
#         mem_mb=get_mem_mb,
#     shell:
#         """
#         export LC_CTYPE=en_US.UTF-8
#         export LC_ALL=en_US.UTF-8
#         workflow/scripts/postprocessing/group_nearby_calls_of_same_AF_and_generate_output_table.pl {input.calls}  > {output.calls}
#         """


rule postprocessing_sv_group_table:
    input:
        calls="{folder}/{sample}/mosaiclassifier/postprocessing/merge/{method}.tsv",
    output:
        grouptrack="{folder}/{sample}/mosaiclassifier/postprocessing/group-table/{method}.tsv",
    log:
        "{folder}/log/mosaiclassifier/postprocessing/group-table/{sample}/{method}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        workflow/scripts/postprocessing/create-sv-group-track.py {input.calls}  > {output.grouptrack}
        """


# rule postprocessing_sv_group_table:
#     input:
#         calls="{folder}/mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#     output:
#         grouptrack="{folder}/mosaiclassifier/postprocessing/group-table/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#     log:
#         "{folder}/log/mosaiclassifier/postprocessing/group-table/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.log",
#     conda:
#         "../envs/mc_base.yaml"
#     resources:
#         mem_mb=get_mem_mb,
#     shell:
#         """
#         PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
#         workflow/scripts/postprocessing/create-sv-group-track.py {input.calls}  > {output.grouptrack}
#         """


rule filter_calls:
    input:
        # inputcalls = "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
        inputcalls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filterFALSE.tsv",
        # mergedcalls = "mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
        mergedcalls="{folder}/{sample}/mosaiclassifier/postprocessing/merge/{method}.tsv",
    output:
        # calls = "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterTRUE.tsv"
        calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filterTRUE.tsv",
    log:
        "{folder}/log/mosaiclassifier/sv_calls/{sample}/{method}_filterTRUE.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        workflow/scripts/postprocessing/apply_filter.py {input.inputcalls} {input.mergedcalls} > {output.calls}
        """


# rule filter_calls:
#     input:
#         # inputcalls = "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
#         inputcalls="{folder}/mosaiclassifier/sv_calls/{sample}/{method}_filterFALSE.tsv",
#         # mergedcalls = "mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
#         mergedcalls="{folder}/mosaiclassifier/postprocessing/merge/{sample}/{method}.tsv",
#     output:
#         # calls = "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterTRUE.tsv"
#         calls="{folder}/mosaiclassifier/sv_calls/{sample}/{method}_filterTRUE.tsv",
#     log:
#         "{folder}/log/mosaiclassifier/sv_calls/{sample}/{method}_filterTRUE.log",
#     conda:
#         "../envs/mc_base.yaml"
#     resources:
#         mem_mb=get_mem_mb,
#     shell:
#         """
#         PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
#         workflow/scripts/postprocessing/apply_filter.py {input.inputcalls} {input.mergedcalls} > {output.calls}
#         """
