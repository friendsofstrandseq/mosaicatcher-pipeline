
################################################################################
# PostProcessing                                                               #
################################################################################


# DOCME : perl in conda

# CHECKME : segdup file only in hg38
rule postprocessing_filter:
    input: 
        calls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv"
    output: 
        calls = config["output_location"] + "mosaiclassifier/postprocessing/filter/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv"
    conda: 
        "../envs/mc_base.yaml"
    shell:
        """
        export LC_CTYPE=en_US.UTF-8 
        export LC_ALL=en_US.UTF-8 
        scripts/postprocessing/filter_MosaiCatcher_calls.pl {input.calls} {config[segdups]} > {output.calls}
        """

rule postprocessing_merge:
    input: 
        calls = config["output_location"] + "mosaiclassifier/postprocessing/filter/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv"
    output: 
        calls = config["output_location"] + "mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv"
    conda: 
        "../envs/mc_base.yaml"
    shell:
        """
        export LC_CTYPE=en_US.UTF-8 
        export LC_ALL=en_US.UTF-8 
        scripts/postprocessing/group_nearby_calls_of_same_AF_and_generate_output_table.pl {input.calls}  > {output.calls}
        """


rule postprocessing_sv_group_table:
    input: 
        calls = config["output_location"] + "mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv"
    output: 
        grouptrack = config["output_location"] + "mosaiclassifier/postprocessing/group-table/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv"
    conda: 
        "../envs/mc_base.yaml"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        scripts/postprocessing/create-sv-group-track.py {input.calls}  > {output.grouptrack}
        """



rule filter_calls:
    input: 
        # inputcalls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterFALSE.tsv",
        inputcalls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}_filterFALSE.tsv",
        # mergedcalls = config["output_location"] + "mosaiclassifier/postprocessing/merge/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}.tsv",
        mergedcalls = config["output_location"] + "mosaiclassifier/postprocessing/merge/{sample}/{method}.tsv",
    output: 
        # calls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/simpleCalls_llr{llr}_poppriors{pop_priors}_haplotags{use_haplotags}_gtcutoff{gtcutoff}_regfactor{regfactor}_filterTRUE.tsv"
        calls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}_filterTRUE.tsv"
    conda: 
        "../envs/mc_base.yaml"
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        scripts/postprocessing/apply_filter.py {input.inputcalls} {input.mergedcalls} > {output.calls}
        """

