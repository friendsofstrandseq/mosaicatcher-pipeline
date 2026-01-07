if config["multistep_normalisation"] is True and config["window"] == 200000:

    rule ashleys_library_size_normalisation:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
            info_raw="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            counts_scaled="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.gz",
        log:
            "{folder}/{sample}/log/counts_scaling/{sample}.log",
        params:
            gc_min_reads=config["multistep_normalisation_options"]["min_reads_cell"],
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../../envs/rtools.yaml"
        script:
            "../../scripts/ashleys/GC/library_size_normalisation.R"

    rule ashleys_GC_correction:
        input:
            counts_scaled="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.gz",
        output:
            counts_scaled_gc="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.gz",
            plot=report(
                "{folder}/{sample}/plots/multistep_normalisation/GC_correction_lowess.png",
                category="GC analysis",
                subcategory="{sample}",
                labels={
                    "Sample": "{sample}",
                    "Plot Type": "Lowess",
                },
            ),
        log:
            "{folder}/{sample}/log/multistep_normalisation/{sample}.log",
        params:
            gc_matrix=ancient(
                "workflow/data/GC/{assembly}.GC_matrix.txt.gz".format(
                    assembly=config["reference"]
                )
            ),
            # gc_matrix=ancient("workflow/data/GC/GC_matrix_200000.txt"),
            gc_min_reads=config["multistep_normalisation_options"]["min_reads_bin"],
            gc_n_subsample=config["multistep_normalisation_options"]["n_subsample"],
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../../envs/rtools.yaml"
        script:
            "../../scripts/ashleys/GC/GC_correction.R"

    rule ashleys_VST_correction:
        input:
            counts_scaled_gc="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.gz",
        output:
            counts_scaled_gc_vst="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.gz",
            plot=report(
                "{folder}/{sample}/plots/multistep_normalisation/GC_correction_VST_hist.png",
                category="GC analysis",
                subcategory="{sample}",
                labels={
                    "Sample": "{sample}",
                    "Plot Type": "VST",
                },
            ),
        log:
            "{folder}/{sample}/log/VST_correction/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../../envs/rtools.yaml"
        script:
            "../../scripts/ashleys/GC/variance_stabilizing_transformation.R"

    rule ashleys_reformat_ms_norm:
        input:
            "{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.gz",
        output:
            "{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.reformat.gz",
        log:
            "{folder}/{sample}/log/reformat_ms_norm/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../../scripts/ashleys/utils/reformat_ms_norm.py"

    rule ashleys_populate_counts_GC:
        input:
            bin_bed=ancient(select_binbed),
            counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.reformat.gz",
        output:
            populated_counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.populated.gz",
        log:
            "{folder}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../../scripts/ashleys/utils/populated_counts_for_qc_plot.py"

    rule ashleys_plot_mosaic_gc_norm_counts:
        input:
            counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.populated.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            # "{folder}/{sample}/plots/counts/CountComplete.normalised.pdf",
            report(
                "{folder}/{sample}/plots/counts/CountComplete.normalised.pdf",
                category="Mosaic Counts",
                subcategory="{sample}",
                labels={"Cell": "ALL", "Type": "normalised"},
            ),
        log:
            "{folder}/{sample}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../../envs/rtools.yaml"
        params:
            mouse_assembly=True if config["reference"] == "mm10" else False,
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """
