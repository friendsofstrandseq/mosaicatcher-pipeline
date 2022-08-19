################################################################################
# Ploidy estimation                                                            #
################################################################################

# if int(config["window"]) in [50000, 100000, 200000]:


rule estimate_ploidy:
    input:
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
    output:
        "{output_folder}/ploidy/{sample}/ploidy_detailled.txt",
    log:
        "{output_folder}/log/ploidy/{sample}/ploidy_detailled.log",
    threads: 48
    resources:
        mem_mb=get_mem_mb,
    params:
        # TODO move this to config
        merge_window=1000000,
        shift_step=500000,
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
            --output {output}
        """


checkpoint summarise_ploidy:
    input:
        ploidy="{output_folder}/ploidy/{sample}/ploidy_detailled.txt",
    output:
        summary="{output_folder}/ploidy/{sample}/ploidy_summary.txt",
    log:
        "{output_folder}/ploidy/{sample}/ploidy_summary.log",
    run:
        df = pd.read_csv(input.ploidy, sep="\t").groupby("#chrom")[
            "ploidy_estimate"
        ].describe()
        df.to_csv(output.summary, sep="\t")
        # haploid_chroms = df.loc[df["50%"] == 1, "#chrom"].values.tolist()
        # if not haploid_chroms:
        #     haploid_chroms = "None"
        # log[0].write("")




rule ploidy_bcftools:
    input:
        "{output_folder}/ploidy/{sample}/ploidy_detailled.txt",
    output:
        "{output_folder}/ploidy/{sample}/ploidy_bcftools.txt",
    log:
        "{output_folder}/ploidy/{sample}/ploidy_bcftools.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/ploidy_bcftools.py"
