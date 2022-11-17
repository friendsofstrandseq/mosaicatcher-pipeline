
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
        tsv="{folder}/{sample}/stats/stats-merged.tsv",
    log:
        tsv="{folder}/log/stats/{sample}/stats-merged.tsv",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "(head -n1 {input.tsv[0]} && (tail -n1 -q {input.tsv} | sort -k1) ) > {output}"


rule transpose_table:
    input:
        "{folder}/{sample}/stats/stats-merged.tsv",
    output:
        html=report(
            "{folder}/{sample}/stats/stats-merged.html",
            category="Stats",
            subcategory="{sample}",
            labels={
                "Sample": "{sample}",
            }
        )   

    run:
        import pandas as pd
        import os
        df = pd.read_csv(input[0], sep="\t")
        df["callset"] = df["callset"].apply(lambda r: os.path.basename(r).replace(".tsv", ""))
        df = df.set_index("callset")
        df = df.fillna(0).T.reset_index()
        pd.options.display.float_format = '{:,.1f}'.format
        df_out = write_to_html_file(df, wildcards.sample)
        with open(output.html, 'w') as o:
            o.write(df_out)