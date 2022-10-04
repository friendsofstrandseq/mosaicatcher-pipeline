# from workflow.scripts.utils.utils import get_mem_mb

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"


envvars:
    "OPENBLAS_NUM_THREADS",


################################################################################
# Plots                                                                        #
################################################################################

# Load rules only if plot is enabled [True] in config file
# if config["plot"] is True:


rule plot_mosaic_counts:
    """
    rule fct: Plot function of read counts for each bam file
    input: mosaic count outputs (counts & info)
    output: Generate figure based on couting results
    """
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        info="{output_folder}/counts/{sample}/{sample}.info",
    output:
        "{output_folder}/plots/{sample}/counts/CountComplete.classic.pdf",
    log:
        "{output_folder}/log/plot_mosaic_counts/{sample}.log",
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
        expand("{{output_folder}}/plots/{{sample}}/counts/CountComplete.{plottype}.pdf", plottype=plottype_counts),
    output:
        report(
            "{output_folder}/plots/{sample}/counts_{plottype}/{cell}.{i, \d+}.pdf",
            caption="../report/mosaic_counts.rst",
            category="Mosaic counts",
            subcategory="{sample}",
            labels={"Cell": "{cell}", "Nb": "{i}", "Type": "{plottype}"},
        ),
    log:
        "{output_folder}/log/plots/{sample}/counts_{plottype}/{cell}.{i, \d+}.log",
    conda:
        "../envs/mc_base.yaml"
    # params:
    #     config_df="{output_folder}/config/config_df.tsv".format(output_folder=config["output_location"]),
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/dividing_pdf.py"


rule tmp_merge_divide:
    input:
        aggregate_cells_count_plot,
    output:
        touch("{output_folder}/plots/{sample}/divide_plots/{sample}.txt"),
    log:
        "{output_folder}/log/final_blank_results/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "touch {output}"


rule final_results:
    input:
        get_all_plots,
    output:
        "{output_folder}/plots/{sample}/final_results/{sample}.txt",
    log:
        "{output_folder}/log/final_blank_results/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "touch {output}"


rule plot_SV_consistency_barplot:
    input:
        sv_calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
    output:
        barplot_bypos=report(
            "{output_folder}/plots/{sample}/sv_consistency/{method}_filter{filter}.consistency-barplot-bypos.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            labels={
                "Barplot type": "By position",
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
        barplot_byaf=report(
            "{output_folder}/plots/{sample}/sv_consistency/{method}_filter{filter}.consistency-barplot-byaf.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            labels={
                "Barplot type": "By AF",
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
    log:
        "{output_folder}/log/plot_SV_consistency/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/sv_consistency_barplot.snakemake.R"


rule plot_clustering:
    input:
        sv_calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
        binbed="workflow/data/bin_200kb_all.bed",
    output:
        position=report(
            "{output_folder}/plots/{sample}/sv_clustering/{method}-filter{filter}-position.pdf",
            category="SV Clustering",
            subcategory="{sample}",
            labels={
                "method": "{method}",
                "filter": "{filter}",
            },
        ),
    log:
        "{output_folder}/log/plot_clustering/{sample}/{method}_filter{filter}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/plot-clustering.snakemake.R"


rule plot_SV_calls:
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
        calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
        complex_calls="{output_folder}/mosaiclassifier/complex/{sample}/{method}_filter{filter}.tsv",
        strand="{output_folder}/strandphaser/{sample}/StrandPhaseR_final_output.txt",
        segments="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
        scsegments="{output_folder}/segmentation/{sample}/Selection_singleseg.txt",
        grouptrack="{output_folder}/mosaiclassifier/postprocessing/group-table/{sample}/{method}.tsv",
    output:
        report(
            "{output_folder}/plots/{sample}/sv_calls/{method}_filter{filter}/{chrom}.pdf",
            category="SV Calls",
            subcategory="{sample}",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chrom": "{chrom}",
            },
        ),
    log:
        "{output_folder}/log/plot_SV_calls/{sample}/{method}_filter{filter}/{chrom}.log",
    conda:
        # "../envs/dev/plot_sv_calls.yaml"
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


rule plot_ploidy:
    input:
        ploidy_detailled="{output_folder}/ploidy/{sample}/ploidy_detailled.txt",
    output:
        report(
            "{output_folder}/plots/{sample}/ploidy/{sample}.pdf",
            category="Ploidy",
            labels={"Sample": "{sample}"}
        ),
    log:
        "{output_folder}/log/plot_ploidy/{sample}.log",
    conda:
        "../envs/python_plots.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/ploidy_plot.py"

if config["GC_analysis"] is True:

    rule plot_mosaic_gc_norm_counts:
    """
    rule fct: Plot function of read counts for each bam file
    input: mosaic count outputs (counts & info)
    output: Generate figure based on couting results
    """
    input:
        counts="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.GC.scaled.txt.gz",
        info="{output_folder}/counts/{sample}/{sample}.info",
    output:
        "{output_folder}/plots/{sample}/counts/CountComplete.GC_corrected.pdf",
    log:
        "{output_folder}/log/plot_mosaic_counts/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
        """

    rule alfred_plot:
        input:
            table = "{output_folder}/alfred/{sample}.table"
        output:
            gcdist_plot = report(
                "{output_folder}/plots/{sample}/alfred/gc_dist.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Type": "GC distribution"}
            ),
            gcdevi_plot = report(
                "{output_folder}/plots/{sample}/alfred/gc_devi.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Type": "GC deviation"}
            ),
            # gcdist_plot = "{output_folder}/plots/{sample}/alfred/gc_dist.png",
            # gcdevi_plot = "{output_folder}/plots/{sample}/alfred/gc_devi.png",
        log:
            "{output_folder}/log/alfred_plot/{sample}.log"
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/plotting/gc.R"


# rule datavzrd_variants_calls:
#     input:
#         config = "workflow/datavzrd/datavzrd_sv_calls.yaml",
#     output:
#         report(
#             directory("results/datavzrd-report"),
#             htmlindex="index.html",
#         ),
#     conda:
#         "datavzrd"
#     shell:
#         "datavzrd {input.config} --output {output}"


# rule plot_SV_calls:
#     input:
#         counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
#         calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
#         complex_calls="{output_folder}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.complex.tsv",
#         strand="{output_folder}/strandphaser/{sample}/StrandPhaseR_final_output.txt",
#         segments="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
#         scsegments="{output_folder}/segmentation/{sample}/Selection_singleseg.txt",
#         grouptrack="{output_folder}/mosaiclassifier/postprocessing/group-table/{sample}/{method}.tsv",
#     output:
#         report(
#             "{output_folder}/plots/{sample}/sv_calls/{method}_filter{filter}.{chrom}.pdf",
#             category="SV Calls",
#             subcategory="{sample}",
#             labels={
#                 "method": "{method}",
#                 "filter": "{filter}",
#                 "Chrom": "{chrom}",
#             },
#         ),
#     log:
#         "{output_folder}/log/plot_SV_calls/{sample}/{method}_filter{filter}.{chrom}.log",
#     conda:
#         "../envs/rtools.yaml"
#     resources:
#         mem_mb=get_mem_mb,
#     shell:
#         """
#         Rscript workflow/scripts/plotting/plot-sv-calls.R \
#             segments={input.segments} \
#             singlecellsegments={input.scsegments} \
#             strand={input.strand} \
#             complex={input.complex_calls} \
#             groups={input.grouptrack} \
#             calls={input.calls} \
#             {input.counts} \
#             {wildcards.chrom} \
#             {output} > {log} 2>&1
#         """
