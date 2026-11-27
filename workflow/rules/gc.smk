if config["multistep_normalisation"] is True and config["window"] == 200000:

    # Implementation of the library size, GC and VST steps:
    #   strandtools: Python, strandtools package (git.embl.org/cosenza/strandtools)
    #   legacy: the former R scripts. DEPRECATED, kept only to compare outputs:
    #     GC_correction.R amplifies the GC bias instead of correcting it (#140).
    # Both read and write the same files: the rules downstream do not change.
    ms_norm_method = config["multistep_normalisation_options"].get(
        "method", "strandtools"
    )
    assert ms_norm_method in (
        "strandtools",
        "legacy",
    ), f"multistep_normalisation_options.method must be 'strandtools' or 'legacy', got '{ms_norm_method}'"
    if ms_norm_method == "legacy":
        from snakemake.logging import logger

        logger.warning(
            "WARNING: multistep_normalisation_options.method=legacy is deprecated. "
            "The R GC correction (GC_correction.R) amplifies the GC bias instead of "
            "correcting it (https://github.com/friendsofstrandseq/mosaicatcher-pipeline/issues/140): "
            "normalised counts and plots are biased. Use method=strandtools."
        )
    ms_norm_ext = "py" if ms_norm_method == "strandtools" else "R"
    ms_norm_env = (
        "../envs/strandtools.yaml"
        if ms_norm_method == "strandtools"
        else "../envs/rtools.yaml"
    )

    rule library_size_normalisation:
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
            ms_norm_env
        script:
            f"../scripts/GC/library_size_normalisation.{ms_norm_ext}"

    rule GC_correction:
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
            seed=config["multistep_normalisation_options"].get("seed", 42),
        resources:
            mem_mb=get_mem_mb,
        conda:
            ms_norm_env
        script:
            f"../scripts/GC/GC_correction.{ms_norm_ext}"

    rule VST_correction:
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
        params:
            # grid searched for the dispersion (strandtools only)
            phi_min=0.001,
            phi_max=1.0,
            phi_steps=50,
        resources:
            mem_mb=get_mem_mb,
        conda:
            ms_norm_env
        script:
            f"../scripts/GC/variance_stabilizing_transformation.{ms_norm_ext}"

    rule reformat_ms_norm:
        input:
            "{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.gz",
        output:
            "{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.reformat.gz",
        log:
            "{folder}/log/reformat_ms_norm/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../scripts/utils/reformat_ms_norm.py"

    rule populate_counts_GC:
        input:
            bin_bed=get_bin_bed_file(),
            counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.reformat.gz",
        output:
            populated_counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.populated.gz",
        log:
            "{folder}/log/populate_counts_GC/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../scripts/utils/populated_counts_for_qc_plot.py"

    rule plot_mosaic_gc_norm_counts:
        input:
            counts="{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.populated.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            report(
                "{folder}/{sample}/plots/counts/CountComplete.normalised.pdf",
                category="Mosaic Counts",
                subcategory="{sample}",
                labels={"Cell": "ALL", "Type": "normalised"},
            ),
        log:
            "{folder}/{sample}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """
