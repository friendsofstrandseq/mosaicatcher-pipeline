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

    checkpoint mosaic_count:
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
            checks=lambda wc: expand(
                "{folder}/{sample}/checks/{cell}.sm_check.ok",
                folder=config["data_location"],
                sample=wc.sample,
                cell=bam_per_sample_local[str(wc.sample)],
            ),
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

    rule populate_counts:
        input:
            bin_bed=ancient("workflow/data/bin_200kb_all.bed"),
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        output:
            populated_counts="{folder}/{sample}/counts/{sample}.txt.populated.gz",
        log:
            "{folder}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../scripts/utils/populated_counts_for_qc_plot.py"

    if config["input_bam_legacy"] is True:

        rule selected_cells:
            input:
                path=ancient("{folder}/{sample}"),
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
        select_labels,
    output:
        "{folder}/{sample}/config/labels.tsv",
    log:
        "{folder}/log/copy_labels/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "cp {input} {output}"


rule symlink_selected_bam:
    input:
        bam="{folder}/{sample}/bam/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    output:
        bam="{folder}/{sample}/selected/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/selected/{cell}.sort.mdup.bam.bai",
    log:
        "{folder}/log/symlink_selected_bam/{sample}/{cell}.log",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/symlink_selected_bam.py"


rule remove_unselected_bam:
    input:
        labels="{folder}/{sample}/cell_selection/labels.tsv",
        bam=unselected_input_bam,
        bai=unselected_input_bai,
    output:
        touch("{folder}/{sample}/config/remove_unselected_bam.ok"),
    log:
        "{folder}/{sample}/log/remove_unselected_bam.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        rm {input.bam} {input.bai}
        """


rule remove_unselected_bam_empty:
    output:
        touch("{folder}/{sample}/config/remove_unselected_bam_empty.ok"),
    log:
        "{folder}/{sample}/log/remove_unselected_bam_empty.log",


checkpoint filter_bad_cells_from_mosaic_count:
    input:
        info_raw="{folder}/{sample}/counts/{sample}.info_raw",
        counts_sort=select_counts_for_SV_calling,
        labels="{folder}/{sample}/config/labels.tsv",
    output:
        info="{folder}/{sample}/counts/{sample}.info",
        info_removed="{folder}/{sample}/counts/{sample}.info_rm",
        counts="{folder}/{sample}/counts/{sample}.txt.filter.gz",
    log:
        "{folder}/log/filter_bad_cells_from_mosaic_count/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/utils/filter_bad_cells.py"


if (
    config["hgsvc_based_normalized_counts"] is True
    and (config["window"] in [50000, 100000, 200000])
    and (config["reference"] == "hg38")
):

    rule merge_blacklist_bins_for_norm:
        input:
            norm=ancient("workflow/data/normalization/{reference}/HGSVC.{window}.txt"),
            whitelist=ancient("workflow/data/normalization/inversion-whitelist.tsv"),
        output:
            merged="{folder}/{sample}/normalizations/{reference}/HGSVC.{window}.merged.tsv",
        log:
            "{folder}/log/merge_blacklist_bins/{sample}/{reference}/HGSVC.{window}.merged.tsv",
        params:
            window=config["window"],
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            workflow/scripts/normalization/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size {params.window} --output {output.merged}
            """

else:

    rule merge_blacklist_bins:
        input:
            norm=ancient(
                "workflow/data/arbigent/normalization/{reference}/HGSVC.{window}.txt"
            ),
        output:
            merged="{folder}/{sample}/normalizations/{reference}/HGSVC.{window}.merged.tsv",
        log:
            "{folder}/log/merge_blacklist_bins_arbigent/{sample}/{reference}/HGSVC.{window}.merged.tsv",
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            cp {input.norm} {ouput.merged}
            """


if config["blacklist_regions"] is True:

    rule normalize_counts:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.filter.gz",
            norm=lambda wc: expand(
                "{folder}/{sample}/normalizations/{reference}/HGSVC.{window}.merged.tsv",
                folder=config["data_location"],
                sample=wc.sample,
                reference=config["reference"],
                window=config["window"],
            ),
        output:
            "{folder}/{sample}/counts/{sample}.txt.gz",
        log:
            "{folder}/log/normalize_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        params:
            normalisation_type=config["hgsvc_based_normalized_counts"],
        shell:
            """
            Rscript workflow/scripts/normalization/normalize.R {input.counts} {input.norm} {output} {params.normalisation_type} 2>&1 > {log}
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
    resources:
        mem_mb=get_mem_mb,
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
