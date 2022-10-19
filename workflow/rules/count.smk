if config["ashleys_pipeline"] is False:

    rule generate_exclude_file_for_mosaic_count:
        input:
            bam=lambda wc: expand(
                "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
                folder=config["data_location"],
                sample=wc.sample,
                cell=bam_per_sample_local[str(wc.sample)],
            ),
        output:
            "{folder}/{sample}/config/chroms_to_exclude.txt",
        log:
            "{folder}/log/config_output/{sample}/exclude_file.log",
        params:
            chroms=config["chromosomes"],
        conda:
            "../envs/mc_base.yaml"
        script:
            "../scripts/utils/generate_exclude_file.py"

    rule mosaic_count:
        input:
            bam=lambda wc: expand(
                "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
                folder=config["data_location"],
                sample=wc.sample,
                cell=bam_per_sample_local[str(wc.sample)],
            ),
            bai=lambda wc: expand(
                "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
                folder=config["data_location"],
                sample=wc.sample,
                cell=bam_per_sample_local[str(wc.sample)],
            ),
            excl="{folder}/{sample}/config/chroms_to_exclude.txt",
        output:
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        log:
            "{folder}/log/counts/{sample}/mosaic_count.log",
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


    if config["input_old_behavior"] is True:

        rule selected_cells:
            input:
                path="{folder}/{sample}",
            output:
                "{folder}/{sample}/cell_selection/labels.tsv",
            log:
                "{folder}/log/{sample}/selected_cells/labels.log",
            conda:
                "../envs/mc_base.yaml"
            script:
                "../scripts/utils/handle_input_old_behavior.py"


    else:

        rule touch_labels:
            output:
                "{folder}/{sample}/cell_selection/labels.tsv",
            log:
                "{folder}/log/{sample}/touch_labels/labels.log",
            conda:
                "../envs/mc_base.yaml"
            shell:
                "echo 'cell\tprobability\tprediction' > {output}"


rule copy_labels:
    input:
        lambda wc: expand(
            "{folder}/{sample}/cell_selection/labels.tsv",
            folder=config["data_location"],
            sample=wc.sample,
        ),
    output:
        "{folder}/{sample}/config/labels.tsv",
    log:
        "{folder}/log/copy_labels/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "cp {input} {output}"


checkpoint filter_bad_cells_from_mosaic_count:
    input:
        info_raw="{folder}/{sample}/counts/{sample}.info_raw",
        counts_sort="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        labels="{folder}/{sample}/config/labels.tsv",
    output:
        info="{folder}/{sample}/counts/{sample}.info",
        info_removed="{folder}/{sample}/counts/{sample}.info_rm",
        counts="{folder}/{sample}/counts/{sample}.txt.filter.gz",
    log:
        "{folder}/log/filter_bad_cells_from_mosaic_count/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/filter_bad_cells.py"


if ((config["window"] in [50000, 100000, 200000])
    and (config["reference"] == "hg38")
    and (config["normalized_counts"] is True)):

    rule merge_blacklist_bins:
        input:
            norm="workflow/data/normalization/HGSVC.{window}.txt",
            whitelist="workflow/data/normalization/inversion-whitelist.tsv",
        output:
            merged="{folder}/{sample}/normalizations/HGSVC.{window}.merged.tsv",
        log:
            "{folder}/log/normalizations/{sample}/HGSVC.{window}.merged.tsv"
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            workflow/scripts/normalization/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2>> {log}
            """

    rule normalize_counts:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.filter.gz",
            norm=expand(
                "{folder}/{sample}/normalizations/HGSVC.{window}.merged.tsv",
                folder=config["data_location"],
                sample=samples,
                window=config["window"],
            ),
        output:
            "{folder}/{sample}/counts/{window}.txt.gz",
        log:
            "{folder}/log/normalize_counts/{sample}_{window}.log",
        conda:
            "../envs/rtools.yaml"
        shell:
            """
            Rscript workflow/scripts/normalization/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
            """
else:
    rule cp_mosaic_count:
        input:
            "{folder}/{sample}/counts/{sample}.txt.filter.gz",
        output:
            "{folder}/{sample}/counts/{sample}.txt.gz",
        log:
            "{folder}/log/counts/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        shell:
            "cp {input} {output}"
rule sort_counts:
    input:
        "{folder}/{sample}/counts/{sample}.txt.gz",
    output:
        "{folder}/{sample}/counts/{sample}.txt.sort.gz",
    log:
        "{folder}/log/sort_counts/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/sort_counts.py"
rule extract_single_cell_counts:
    input:
        info="{folder}/{sample}/counts/{sample}.info",
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
    output:
        "{folder}/{sample}/counts/counts-per-cell/{cell}.txt.percell.gz",
    log:
        "{folder}/log/counts/{sample}/counts-per-cell/{cell}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        zcat {input.counts} | awk -v name={wildcards.cell} '(NR==1) || $5 == name' | gzip > {output}
        """
