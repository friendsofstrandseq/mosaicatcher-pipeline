import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"


envvars:
    "OPENBLAS_NUM_THREADS",


if config["ashleys_pipeline"] is False:

    rule plot_mosaic_counts:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.populated.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            # "{folder}/{sample}/plots/counts/CountComplete.raw.pdf",
            report(
                "{folder}/{sample}/plots/counts/CountComplete.raw.pdf",
                category="Mosaic counts",
                subcategory="{sample}",
                labels={"Cell": "ALL", "Type": "raw"},
            ),
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
        "{folder}/{sample}/plots/counts/CountComplete.raw.pdf",
    output:
        report(
            "{folder}/{sample}/plots/counts_raw/{cell}.{i, \d+}.pdf",
            caption="../report/mosaic_counts.rst",
            category="Mosaic counts cellwise",
            subcategory="{sample}",
            labels={"Cell": "{cell}", "Nb": "{i}", "Type": "raw"},
        ),
    log:
        "{folder}/log/{sample}/plots/counts_raw/{cell}.{i, \d+}.log",
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
        sv_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv"
        ),
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
        sv_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv"
        ),
        binbed=ancient("workflow/data/bin_200kb_all.bed"),
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
        sv_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv"
        ),
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
        sv_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/{method}_filter{filter}.tsv"
        ),
        binbed=ancient(select_binbed),
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
        complex_calls=(
            "{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv"
        ),
        strand="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
        segments="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        scsegments="{folder}/{sample}/segmentation/Selection_singleseg.txt",
        grouptrack=(
            "{folder}/{sample}/mosaiclassifier/postprocessing/group-table/{method}.tsv"
        ),
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
        complex_calls=(
            "{folder}/{sample}/mosaiclassifier/complex/{method}_filter{filter}.tsv"
        ),
        strand="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
        segments="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        scsegments="{folder}/{sample}/segmentation/Selection_singleseg.txt",
        grouptrack=(
            "{folder}/{sample}/mosaiclassifier/postprocessing/group-table/{method}.tsv"
        ),
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


rule scTRIP_multiplot:
    input:
        install_check="workflow/config/scTRIP_multiplot.ok",
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        haplotag_bam="{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
        sv_counts="{folder}/{sample}/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv",
    output:
        figure=report(
            "{folder}/{sample}/plots/scTRIP_multiplot/{cell}/{chrom}.pdf",
            category="scTRIP multiplot",
            subcategory="{sample}",
            labels={"Cell": "{cell}", "Chrom": "{chrom}"},
        ),
    log:
        "{folder}/log/scTRIP_multiplot/{sample}/{cell}/{chrom}.log",
    conda:
        "../envs/rtools.yaml"
    container:
        None
    resources:
        mem_mb=get_mem_mb,
    shell:
        "LC_CTYPE=C Rscript workflow/scripts/plotting/scTRIP_multiplot/scTRIP_multiplot_run.R {input.counts} {input.haplotag_bam} {input.sv_counts} {wildcards.chrom} {wildcards.cell} {output.figure} > {log} 2>&1"


rule scTRIP_multiplot_aggr:
    input:
        aggregate_cells_scTRIP_multiplot,
    output:
        touch("{folder}/{sample}/plots/scTRIP_multiplot_aggr.ok"),
    log:
        "{folder}/log/scTRIP_multiplot_aggr/{sample}.log",
    resources:
        mem_mb=get_mem_mb,


rule ucsc_genome_browser_file:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        stringent_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv"
        ),
        lenient_calls=(
            "{folder}/{sample}/mosaiclassifier/sv_calls/lenient_filterFALSE.tsv"
        ),
    output:
        "{folder}/{sample}/plots/UCSC/{sample}.bedUCSC.gz",
    log:
        "{folder}/log/ucsc_genome_browser_file/{sample}.bedUCSC.gz",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "python workflow/scripts/plotting/ucsc_vizu.py {input.counts} {input.stringent_calls} {input.lenient_calls} {output} > {log}"


rule split_ucsc_into_individual_tracks:
    input:
        ucsc_file="{folder}/{sample}/plots/UCSC/{sample}.bedUCSC.gz",
    output:
        output_dir=directory("{folder}/{sample}/plots/IGV/SPLITTED"),
    log:
        "{folder}/log/split_ucsc_into_individual_tracks/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "sh workflow/scripts/plotting/split_ucsc_file.sh {input.ucsc_file} {output.output_dir}"


rule generate_igv_session:
    input:
        splitted_files_dir="{folder}/{sample}/plots/IGV/SPLITTED",
    output:
        xml_session="{folder}/{sample}/plots/IGV/{sample}_IGV_session.xml",
    log:
        "{folder}/log/generate_igv_session/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "sh workflow/scripts/plotting/generate_IGV_session.sh {input.splitted_files_dir} {output.xml_session}"
