
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
        tsv=lambda wc: [
            sub_e
            for e in [
                expand(
                    "{folder}/{sample}/stats/{method}_filter{filter}.tsv",
                    folder=config["data_location"],
                    sample=wc.sample,
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
    log:
        tsv="{folder}/log/{sample}/transpose_table.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/stats/transpose_table.py"
