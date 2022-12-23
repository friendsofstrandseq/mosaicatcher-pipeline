rule estimate_ploidy:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.sort.gz",
    output:
        "{folder}/{sample}/ploidy/ploidy_detailled.txt",
    log:
        "{folder}/log/estimate_ploidy/{sample}.log",
    threads: 48
    resources:
        mem_mb=get_mem_mb,
    params:
        # TODO move this to config
        merge_window=1000000,
        shift_step=int(config["window"]) * 5,
        boundary_alpha=0.05,
        max_ploidy=6,
        add_bg_component="to_be_done",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        python workflow/scripts/utils/ploidy_estimator.py --debug \
            --merge-bins-to {params.merge_window} \
            --shift-window-by {params.shift_step} \
            --max-ploidy {params.max_ploidy} \
            --boundary-alpha {params.boundary_alpha} \
            --jobs {threads} \
            --input {input.counts} \
            --output {output} \
            --log {log}
        """


checkpoint summarise_ploidy:
    input:
        ploidy="{folder}/{sample}/ploidy/ploidy_detailled.txt",
    output:
        summary="{folder}/{sample}/ploidy/ploidy_summary.txt",
    log:
        "{folder}/log/ploidy/{sample}/ploidy_summary.log",
    run:
        df = (
            pd.read_csv(input.ploidy, sep="\t")
            .groupby("#chrom")["ploidy_estimate"]
            .describe()
        )
        df.to_csv(output.summary, sep="\t")


rule ploidy_bcftools:
    input:
        "{folder}/{sample}/ploidy/ploidy_detailled.txt",
    output:
        "{folder}/{sample}/ploidy/ploidy_bcftools.txt",
    log:
        "{folder}/log/ploidy/{sample}/ploidy_bcftools.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/ploidy_bcftools.py"
