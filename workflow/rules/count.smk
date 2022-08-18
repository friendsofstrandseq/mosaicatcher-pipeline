rule mosaic_count:
    input:
        bam=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam",
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)],
        ),
        bai=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)],
        ),
        excl="{output_folder}/config_output/{sample}/exclude_file",
    output:
        counts="{output_folder}/counts/{sample}/{sample}.txt.raw.gz",
        info="{output_folder}/counts/{sample}/{sample}.info_raw",
    log:
        "{output_folder}/log/counts/{sample}/mosaic_count.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    params:
        window=config["window"],
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {params.window} \
            {input.bam} \
        > {log} 2>&1
        """


if config["ashleys_pipeline"] is False:

    if config["input_old_behavior"] is True:

        rule selected_cells:
            input:
                path="{input_folder}/{sample}",
            output:
                "{input_folder}/{sample}/cell_selection/labels.tsv",
            log:
                "{input_folder}/log/{sample}/selected_cells/labels.log",
            conda:
                "../envs/mc_base.yaml"
            script:
                "../scripts/utils/handle_input_old_behavior.py"


    else:

        rule touch_labels:
            output:
                "{input_folder}/{sample}/cell_selection/labels.tsv",
            log:
                "{input_folder}/log/{sample}/touch_labels/labels.log",
            conda:
                "../envs/mc_base.yaml"
            shell:
                "echo 'cell\tprobability\tprediction' > {output}"


rule copy_labels:
    input:
        expand(
            "{input_folder}/{sample}/cell_selection/labels.tsv",
            input_folder=config["input_bam_location"],
            sample=samples,
        ),
    output:
        "{output_folder}/cell_selection/{sample}/labels.tsv",
    log:
        "{output_folder}/log/copy_labels/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "cp {input} {output}"


rule order_mosaic_count_output:
    input:
        raw_count="{output_folder}/counts/{sample}/{sample}.txt.raw.gz",
        labels="{output_folder}/cell_selection/{sample}/labels.tsv",
    output:
        "{output_folder}/counts/{sample}/{sample}.txt.sort.gz",
    log:
        "{output_folder}/log/counts/{sample}/{sample}.log",
    run:
        df = pd.read_csv(input.raw_count, compression="gzip", sep="\t")
        df = df.sort_values(by=["sample", "cell", "chrom", "start"])
        df.to_csv(output[0], index=False, compression="gzip", sep="\t")


checkpoint filter_bad_cells_from_mosaic_count:
    input:
        info_raw="{output_folder}/counts/{sample}/{sample}.info_raw",
        counts_sort="{output_folder}/counts/{sample}/{sample}.txt.sort.gz",
        labels="{output_folder}/cell_selection/{sample}/labels.tsv",
    output:
        info="{output_folder}/counts/{sample}/{sample}.info",
        info_removed="{output_folder}/counts/{sample}/{sample}.info_rm",
        counts="{output_folder}/counts/{sample}/{sample}.txt.filter.gz",
    log:
        "{output_folder}/log/filter_bad_cells_from_mosaic_count/{sample}/{sample}.log",
    script:
        "../scripts/utils/filter_bad_cells.py"


if (
    (config["window"] in [50000, 100000, 200000])
    and (config["reference"] == "hg38")
    and (config["normalized_counts"] is True)
):

    rule merge_blacklist_bins:
        input:
            norm="workflow/data/normalization/HGSVC.{{window}}.txt".format(
                config["window"]
            ),
            whitelist="workflow/data/normalization/inversion-whitelist.tsv",
        output:
            merged="{{output_folder}}/normalizations/HGSVC.{{window}}.merged.tsv".format(
                config["output_location"], config["window"]
            ),
        log:
            "{{output_folder}}/log/merge_blacklist_bins/{{window}}.log".format(
                config["output_location"], config["window"]
            ),
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            workflow/scripts/normalization/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2>> {log}
            """

    rule normalize_counts:
        input:
            counts="{output_folder}/counts/{sample}/{sample}.txt.filter.gz",
            norm=expand(
                "{output_folder}/normalizations/HGSVC.{window}.merged.tsv",
                output_folder=config["output_location"],
                window=config["window"],
            ),
        output:
            "{output_folder}/counts/{sample}/{sample}.txt.norm.gz",
        log:
            "{output_folder}/log/normalize_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        shell:
            """
            Rscript workflow/scripts/normalization/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
            """

    rule sort_norm_counts:
        input:
            "{output_folder}/counts/{sample}/{sample}.txt.norm.gz",
        output:
            "{output_folder}/counts/{sample}/{sample}.txt.gz",
        log:
            "{output_folder}/sort_norm_counts/{sample}.log",
        run:
            df = pd.read_csv(input[0], sep="\t", compression="gzip")
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)
            chroms = ["chr{}".format(str(c)) for c in list(range(1, 23))] + ["chrX"]
            df["chrom"] = pd.Categorical(df["chrom"], categories=chroms, ordered=True)
            df.sort_values(by=["cell", "chrom", "start", "end"]).to_csv(
                output[0], compression="gzip", sep="\t", index=False
            )


else:

    rule cp_mosaic_count:
        input:
            "{output_folder}/counts/{sample}/{sample}.txt.filter.gz",
        output:
            "{output_folder}/counts/{sample}/{sample}.txt.gz",
        log:
            "{output_folder}/log/counts/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        shell:
            "cp {input} {output}"


rule extract_single_cell_counts:
    input:
        info="{output_folder}/counts/{sample}/{sample}.info",
        counts="{output_folder}/counts/{sample}/{sample}.txt.gz",
    output:
        "{output_folder}/counts/{sample}/counts-per-cell/{cell}.txt.gz",
    log:
        "{output_folder}/log/counts/{sample}/counts-per-cell/{cell}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        zcat {input.counts} | awk -v name={wildcards.cell} '(NR==1) || $5 == name' | gzip > {output}
        """
