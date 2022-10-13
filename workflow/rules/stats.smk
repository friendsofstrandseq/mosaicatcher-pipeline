# import pandas as pd
# config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# samples = sorted(config_df.Sample.unique().tolist())

# samples = sorted(df_config_files.Sample.unique().tolist())

################################################################################
# Summary statistics on sv calls                                               #
################################################################################


rule summary_statistics:
    input:
        segmentation="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        strandstates="{folder}/{sample}/segmentation/Selection_initial_strand_state",
        sv_calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        complex="{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv",
        merged="{folder}/{sample}/mosaiclassifier/postprocessing/merge/{method}.tsv",
    output:
        tsv="{folder}/{sample}/stats/{method}_filter{filter}.tsv",
    log:
        "{folder}/log/summary_statistics/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/stats/summary_stats.py"


rule aggregate_summary_statistics:
    input:
        tsv=[
            sub_e
            for e in [
                expand(
                    "{folder}/{sample}/stats/{method}_filter{filter}.tsv",
                    folder=config["data_location"],
                    sample=samples,
                    method=method,
                    filter=config["methods"][method]["filter"],
                )
                for method in config["methods"]
            ]
            for sub_e in e
        ],
    output:
        # tsv=report("stats/{sample}/stats-merged.tsv", category="Stats", labels={"Type" : "Complete stats"})
        tsv="{folder}/{sample}/stats/stats-merged.tsv",
    log:
        tsv="{folder}/log/stats/{sample}/stats-merged.tsv",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "(head -n1 {input.tsv[0]} && (tail -n1 -q {input.tsv} | sort -k1) ) > {output}"


# rule summary_statistics_to_pdf:
#     input:
#         "stats/{sample}/stats-merged.tsv"
#     output:
#         report("stats/{sample}/stats-merged.pdf", category="Stats", labels={"Type" : "Complete stats"})
#     run:
#         import pandas as pd
#         df = pd.read_csv(input[0], sep='\t')
#         df.to_pdf
