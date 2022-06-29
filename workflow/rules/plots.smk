# from workflow.scripts.utils.utils import get_mem_mb

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"


envvars:
    "OPENBLAS_NUM_THREADS",


# TODO : option for segmentation : enable/disable back colors

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
        counts="{output}/counts/{sample}/{sample}.txt.gz",
        info="{output}/counts/{sample}/{sample}.info",
    output:
        "{output}/plots/{sample}/counts/CountComplete.pdf",
    log:
        "{output}/log/plot_mosaic_counts/{sample}.log",
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
        "{output}/plots/{sample}/counts/CountComplete.pdf",
    output:
        report(
            "{output}/plots/{sample}/counts/{cell}.{i, \d+}.pdf",
            caption="../report/mosaic_counts.rst",
            category="Mosaic counts",
            subcategory="{sample}",
            labels={"Cell": "{cell}", "Nb": "{i}"},
        ),
    log:
        "{output}/log/plots/{sample}/counts/{cell}.{i, \d+}.log",
    conda:
        "../envs/mc_base.yaml"
    params:
        config_df=config["samples"],
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/dividing_pdf.py"

rule tmp_merge_divide:
    input:
        get_indiv_plots_count(),
    output:
        touch("{output}/plots/{sample}/final_results/{sample}.txt"),
    log:
        "{output}/log/final_blank_results/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "echo {input} > {output}"

rule plot_SV_consistency_barplot:
    input:
        sv_calls="{output}/mosaiclassifier/sv_calls/{sample}/{method}.tsv",
    output:
        barplot_bypos=report(
            "{output}/plots/{sample}/sv_consistency/{method}.consistency-barplot-bypos.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            labels={
                "Barplot type": "By position",
                "method": "{method}",
            },
        ),
        barplot_byaf=report(
            "{output}/plots/{sample}/sv_consistency/{method}.consistency-barplot-byaf.pdf",
            category="SV Consistency",
            subcategory="{sample}",
            labels={
                "Barplot type": "By AF",
                "method": "{method}",
            },
        ),
    log:
        "{output}/log/plot_SV_consistency/{sample}/{method}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/sv_consistency_barplot.snakemake.R"

rule plot_clustering:
    input:
        sv_calls="{output}/mosaiclassifier/sv_calls/{sample}/{method}.tsv",
        binbed="workflow/data/bin_200kb_all.bed",
    output:
        position=report(
            "{output}/plots/{sample}/sv_clustering/{method}-position.pdf",
            category="SV Clustering",
            subcategory="{sample}",
            labels={
                "method": "{method}",
            },
        ),
        chromosome="{output}/plots/{sample}/sv_clustering/{method}-chromosome.pdf",
        # chromosome = report("plots/{sample}/sv_clustering/{method}-chromosome.pdf", category="SV clustering"),
    log:
        "{output}/log/plot_clustering/{sample}/{method}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/plotting/plot-clustering.snakemake.R"
        # binbed = "data/bin_200kb_all.bed",


#        Rscript scripts/plotting/plot-clustering.snakemake.R {input.sv_calls} {input.binbed} {output.position} {output.chromosome}

rule plot_SV_calls:
    input:
        counts="{output}/counts/{sample}/{sample}.txt.gz",
        calls="{output}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
        complex_calls="{output}/mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.complex.tsv",
        strand="{output}/strandphaser/{sample}/StrandPhaseR_final_output.txt",
        segments="{output}/segmentation/{sample}/Selection_jointseg.txt",
        scsegments="{output}/segmentation/{sample}/Selection_singleseg.txt",
        grouptrack="{output}/mosaiclassifier/postprocessing/group-table/{sample}/{method}.tsv",
    output:
        # "plots/{sample}/sv_calls/{method}_filter{filter}.{chrom}.pdf"
        report(
            "{output}/plots/{sample}/sv_calls/{method}_filter{filter}.{chrom}.pdf",
            category="SV Calls",
            subcategory="{sample}",
            labels={
                "method": "{method}",
                "filter": "{filter}",
                "Chrom": "{chrom}",
            },
        ),
    log:
        "{output}/log/plot_SV_calls/{sample}/{method}_filter{filter}.{chrom}.log",
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
        


# rule generate_halo_json:
#     input:
#         counts = "counts/{sample}/{windows}.txt.gz",
#     output:
#         json = "halo/{sample}/{windows}.json.gz",
#     log:
#         "log/generate_halo_json/{sample}/{windows}.{windows}.log"
#     shell:
#         """
#         PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
#         (./utils/counts_to_json.py {input.counts} | gzip > {output.json})
#         """
