import pandas as pd
config_df = pd.read_csv(config["output_location"] + "config/config_df.tsv", sep="\t")
samples = sorted(config_df.Sample.unique().tolist())

################################################################################
# Summary statistics on sv calls                                               #
################################################################################

rule summary_statistics:
    input:
        segmentation = config["output_location"] + 'segmentation/{sample}/Selection_jointseg.txt',
        strandstates = config["output_location"] + 'segmentation/{sample}/Selection_initial_strand_state',
        sv_calls = config["output_location"] + 'mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv',
        complex = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.complex.tsv",
        merged = config["output_location"] + "mosaiclassifier/postprocessing/merge/{sample}/{method}.tsv",
    output:
        tsv = config["output_location"] + 'stats/{sample}/{method}_filter{filter}.tsv',
    log:
        config["output_location"] + 'log/summary_statistics/{sample}/{method}_filter{filter}.log'
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
        shell('scripts/stats/callset_summary_stats.py --segmentation {input.segmentation} --strandstates {input.strandstates} --complex-regions {input.complex} {additional_params} {input.sv_calls}  > {output.tsv} ')

rule aggregate_summary_statistics:
    input:
        tsv=expand(config["output_location"] + "stats/{sample}/{method}.tsv", method=config['methods'], sample=samples),
    output:
        # tsv=report(config["output_location"] + "stats/{sample}/stats-merged.tsv", category="Stats", labels={"Type" : "Complete stats"})
        tsv=config["output_location"] + "stats/{sample}/stats-merged.tsv"
    shell:
        "(head -n1 {input.tsv[0]} && (tail -n1 -q {input.tsv} | sort -k1) ) > {output}"

# rule summary_statistics_to_pdf:
#     input:
#         config["output_location"] + "stats/{sample}/stats-merged.tsv"
#     output:
#         report(config["output_location"] + "stats/{sample}/stats-merged.pdf", category="Stats", labels={"Type" : "Complete stats"})
#     run:
#         import pandas as pd
#         df = pd.read_csv(input[0], sep='\t')
#         df.to_pdf