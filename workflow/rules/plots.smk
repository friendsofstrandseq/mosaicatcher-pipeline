import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"


envvars:
    "OPENBLAS_NUM_THREADS",


if config["ashleys_pipeline"] is False:

    rule plot_mosaic_counts:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            "{folder}/{sample}/plots/counts/CountComplete.classic.pdf",
        log:
            "{folder}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """


rule divide_pdf:
    input:
        "{folder}/{sample}/plots/counts/CountComplete.{plottype}.pdf",
    output:
        report(
            "{folder}/{sample}/plots/counts_{plottype}/{cell}.{i, \d+}.pdf",
            caption="../report/mosaic_counts.rst",
            category="Mosaic counts",
            subcategory="{sample}",
            labels={"Cell": "{cell}", "Nb": "{i}", "Type": "{plottype}"},
        ),
    log:
        "{folder}/log/{sample}/plots/counts_{plottype}/{cell}.{i, \d+}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/dividing_pdf.py"


rule final_results:
    input:
        get_all_plots,
    output:
        "{folder}/{sample}/plots/final_results/{sample}.txt",
    log:
        "{folder}/log/final_blank_results/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "touch {output}"


rule plot_SV_consistency_barplot:
    input:
        sv_calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
    output:
        barplot_bypos=report(
            "{folder}/{sample}/plots/sv_consistency/{method}_filter{filter}.consistency-barplot-bypos.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            caption="../report/sv_consistency.rst",
            labels={
                "Barplot type": "By position",
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
        barplot_byaf=report(
            "{folder}/{sample}/plots/sv_consistency/{method}_filter{filter}.consistency-barplot-byaf.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            caption="../report/sv_consistency.rst",
            labels={
                "Barplot type": "By AF",
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
    log:
        "{folder}/log/plot_SV_consistency/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/sv_consistency_barplot.snakemake.R"


rule plot_clustering:
    input:
        sv_calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        binbed="workflow/data/bin_200kb_all.bed",
    output:
        position=report(
            "{folder}/{sample}/plots/sv_clustering/{method}-filter{filter}-position.pdf",
            category="SV Clustering",
            subcategory="{sample}",
            caption="../report/sv_clustering.rst",
            labels={
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
    log:
        "{folder}/log/plot_clustering/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/sv_heatmap.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/plot-clustering.snakemake.R"


rule plot_clustering_position_dev:
    input:
        sv_calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
    output:
        pdf=report(
            "{folder}/{sample}/plots/sv_clustering_dev/{method}-filter{filter}-position.pdf",
            category="SV Clustering",
            subcategory="{sample}",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chr size scaled": "False",
            },
        ),
        cluster_order_df="{folder}/{sample}/plots/sv_clustering_dev/clustering_{method}-filter{filter}-position.tsv",
    log:
        "{folder}/log/plot_clustering_dev/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/plot_clustering_dev_clean.R"


rule plot_clustering_chromosome_dev:
    input:
        sv_calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        binbed="workflow/data/bin_200kb_all.bed",
        cluster_order_df="{folder}/{sample}/plots/sv_clustering_dev/clustering_{method}-filter{filter}-position.tsv",
    output:
        pdf=report(
            "{folder}/{sample}/plots/sv_clustering_dev/{method}-filter{filter}-chromosome.pdf",
            category="SV Clustering",
            subcategory="{sample}",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chr size scaled": "True",
            },
        ),
    log:
        "{folder}/log/plot_clustering_chromosome_dev/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/plot_clustering_scale_clean.py"


rule plot_SV_calls:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        complex_calls="{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv",
        strand="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
        segments="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        scsegments="{folder}/{sample}/segmentation/Selection_singleseg.txt",
        grouptrack="{folder}/{sample}/mosaiclassifier/postprocessing/group-table/{method}.tsv",
    output:
        report(
            "{folder}/{sample}/plots/sv_calls/{method}_filter{filter}/{chrom}.pdf",
            category="SV Calls",
            subcategory="{sample}",
            caption="../report/sv_calls.rst",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chrom": "{chrom}",
            },
        ),
    log:
        "{folder}/log/plot_SV_calls/{sample}/{method}_filter{filter}/{chrom}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        Rscript workflow/scripts/plotting/plot-sv-calls.R \
            segments={input.segments} \
            singlecellsegments={input.scsegments} \
            strand={input.strand} \
            complex={input.complex_calls} \
            groups={input.grouptrack} \
            calls={input.calls} \
            {input.counts} \
            {wildcards.chrom} \
            {output} > {log} 2>&1
        """


rule plot_SV_calls_dev:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        calls="{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv",
        complex_calls="{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv",
        strand="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
        segments="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        scsegments="{folder}/{sample}/segmentation/Selection_singleseg.txt",
        grouptrack="{folder}/{sample}/mosaiclassifier/postprocessing/group-table/{method}.tsv",
    output:
        report(
            "{folder}/{sample}/plots/sv_calls_dev/{method}_filter{filter}/{chrom}.pdf",
            category="SV Calls",
            subcategory="{sample}",
            caption="../report/sv_calls.rst",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chrom": "{chrom}",
            },
        ),
    log:
        "{folder}/log/plot_SV_calls_dev/{sample}/{method}_filter{filter}/{chrom}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        Rscript workflow/scripts/plotting/plot-sv-calls_dev.R \
            segments={input.segments} \
            singlecellsegments={input.scsegments} \
            strand={input.strand} \
            complex={input.complex_calls} \
            groups={input.grouptrack} \
            calls={input.calls} \
            {input.counts} \
            {wildcards.chrom} \
            {output} > {log} 2>&1
        """


rule plot_ploidy:
    input:
        ploidy_detailled="{folder}/{sample}/ploidy/ploidy_detailled.txt",
    output:
        report(
            "{folder}/{sample}/plots/ploidy/{sample}.pdf",
            category="Ploidy",
            labels={"Sample": "{sample}"},
        ),
    log:
        "{folder}/log/plot_ploidy/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/ploidy_plot.py"
