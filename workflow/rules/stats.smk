# import pandas as pd
# config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# samples = sorted(config_df.Sample.unique().tolist())

# samples = sorted(df_config_files.Sample.unique().tolist())

################################################################################
# Summary statistics on sv calls                                               #
################################################################################


rule summary_statistics:
    input:
        segmentation="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
        strandstates="{output_folder}/segmentation/{sample}/Selection_initial_strand_state",
        sv_calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
        complex="{output_folder}/mosaiclassifier/complex/{sample}/{method}_filter{filter}.tsv",
        merged="{output_folder}/mosaiclassifier/postprocessing/merge/{sample}/{method}.tsv",
    output:
        tsv="{output_folder}/stats/{sample}/{method}_filter{filter}.tsv",
    log:
        "{output_folder}/log/summary_statistics/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/stats/summary_stats.py"


rule aggregate_summary_statistics:
    input:
        tsv=expand(
            "{output_folder}/stats/{sample}/{method}_filter{filter}.tsv",
            output_folder=config["output_location"],
            method=config["methods"],
            sample=samples,
        ),
    output:
        # tsv=report("stats/{sample}/stats-merged.tsv", category="Stats", labels={"Type" : "Complete stats"})
        tsv="{output_folder}/stats/{sample}/stats-merged.tsv",
    log:
        tsv="{output_folder}/log/stats/{sample}/stats-merged.tsv",
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
