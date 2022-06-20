from scripts.utils.utils import get_mem_mb 

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

# TODO : option for segmentation : enable/disable back colors

################################################################################
# Plots                                                                        #
################################################################################

# Load rules only if plot is enabled [True] in config file        
if config["plot"] is True:

    rule plot_mosaic_counts:
        """
        rule fct: Plot function of read counts for each bam file
        input: mosaic count outputs (counts & info)
        output: Generate figure based on couting results
        """
        input:
            counts = config["output_location"] + "counts/{sample}/{sample}.txt.gz",
            info   = config["output_location"] + "counts/{sample}/{sample}.info"
        output:
            config["output_location"] + "plots/{sample}/counts/CountComplete.pdf"
        log:
            config["output_location"] + "log/plot_mosaic_counts/{sample}.log"
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb = get_mem_mb,
        shell:
            """
            Rscript scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """
    # TODO : from shell to script function

    rule divide_pdf:
        input:
            config["output_location"] + "plots/{sample}/counts/CountComplete.pdf"
        output:
            report(
                config["output_location"] + "plots/{sample}/counts/{cell}.{i, \d+}.pdf",
                    caption="../report/mosaic_counts.rst",
                    category="Mosaic counts",
                    subcategory = "{sample}",
                    labels={"Cell" : "{cell}", "Nb" : "{i}"}
            )
        conda:
            "../envs/mc_base.yaml"
        params:
            config_df = config["output_location"] + "config/config_df.tsv"
        resources:
            mem_mb = get_mem_mb,
        script:
            "../scripts/plotting/dividing_pdf.py"


    rule plot_SV_consistency_barplot:
        input:
            sv_calls  = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}.tsv",
        output:
            barplot_bypos = report(config["output_location"] + "plots/{sample}/sv_consistency/{method}.consistency-barplot-bypos.pdf", category="SV Consistency", subcategory="{sample}", labels={"Barplot type" : "By position", "method" : "{method}", }),
            barplot_byaf = report(config["output_location"] + "plots/{sample}/sv_consistency/{method}.consistency-barplot-byaf.pdf", category="SV Consistency", subcategory="{sample}", labels={"Barplot type" : "By AF", "method" : "{method}", }),
        log:
            config["output_location"] + "log/plot_SV_consistency/{sample}/{method}.log"
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb = get_mem_mb,
        script:
            "../scripts/plotting/sv_consistency_barplot.snakemake.R"


    rule plot_clustering:
        input:
            sv_calls  = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}.tsv",
            binbed = "workflow/data/bin_200kb_all.bed",
        output:
            position = report(config["output_location"] + "plots/{sample}/sv_clustering/{method}-position.pdf", category="SV Clustering", subcategory="{sample}", labels={"method" : "{method}", }),
            chromosome = config["output_location"] + "plots/{sample}/sv_clustering/{method}-chromosome.pdf",
            # chromosome = report(config["output_location"] + "plots/{sample}/sv_clustering/{method}-chromosome.pdf", category="SV clustering"),
        log:
            config["output_location"] + "log/plot_clustering/{sample}/{method}.log"
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb = get_mem_mb,
        script:
            "../scripts/plotting/plot-clustering.snakemake.R"
    #        Rscript scripts/plotting/plot-clustering.snakemake.R {input.sv_calls} {input.binbed} {output.position} {output.chromosome}

    rule plot_SV_calls:
        input:
            counts = config["output_location"] + "counts/{sample}/{sample}.txt.gz",
            calls  = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.tsv",
            complex_calls = config["output_location"] + "mosaiclassifier/sv_calls/{sample}/{method}_filter{filter}.complex.tsv",
            strand = config["output_location"] + "strandphaser/{sample}/StrandPhaseR_final_output.txt",
            segments = config["output_location"] + "segmentation/{sample}/Selection_jointseg.txt",
            scsegments = config["output_location"] + "segmentation/{sample}/Selection_singleseg.txt",
            grouptrack = config["output_location"] + "mosaiclassifier/postprocessing/group-table/{sample}/{method}.tsv",
        output:
            # config["output_location"] + "plots/{sample}/sv_calls/{method}_filter{filter}.{chrom}.pdf"
            report(config["output_location"] + "plots/{sample}/sv_calls/{method}_filter{filter}.{chrom}.pdf", category="SV Calls", subcategory="{sample}", labels={"method" : "{method}",  "filter" : "{filter}", "Chrom" : "{chrom}"})
        log:
            config["output_location"] + "log/plot_SV_calls/{sample}/{method}_filter{filter}.{chrom}.log"
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb = get_mem_mb,
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

    # TODO : rule merge_plot_SV_calls ?


    # TODO : from shell to script function



        # rule generate_halo_json:
        #     input:
        #         counts = config["output_location"] + "counts/{sample}/{windows}.txt.gz",
        #     output:
        #         json = config["output_location"] + "halo/{sample}/{windows}.json.gz",
        #     log:
        #         config["output_location"] + "log/generate_halo_json/{sample}/{windows}.{windows}.log"
        #     shell:
        #         """
        #         PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        #         (./utils/counts_to_json.py {input.counts} | gzip > {output.json}) 
        #         """
