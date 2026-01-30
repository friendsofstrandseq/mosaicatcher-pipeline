## Rules dedicated to Strand-Seq count based on mosaic program
## ---------------------------------------------------------------
## generate_exclude_file_for_mosaic_count: generate a list of chromosomes to exclude except the canonical ones
## mosaic_count: mosaic count program to count reads in each bin based on window selected (default: 200kb)
## plot_mosaic_counts: plot QC plots based on counts


rule ashleys_generate_exclude_file_for_mosaic_count:
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
    output:
        excl="{folder}/{sample}/config/chroms_to_exclude.txt",
    log:
        "{folder}/log/config/{sample}/exclude_file.log",
    conda:
        "../../envs/mc_base.yaml"
    params:
        chroms=config["chromosomes"],  # Use config directly for all references
    script:
        "../../scripts/ashleys/utils/generate_exclude_file.py"


checkpoint mosaic_count:
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        bai=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        excl="{folder}/{sample}/config/chroms_to_exclude.txt",
    output:
        counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        info="{folder}/{sample}/counts/{sample}.info_raw",
    log:
        "{folder}/log/counts/{sample}/mosaic_count.log",
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    params:
        window=config["window"],
    resources:
        mem_mb=get_mem_mb_heavy,
        time="24:00:00",
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


rule ashleys_populate_counts:
    input:
        bin_bed=ancient(select_binbed),
        counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
    output:
        populated_counts="{folder}/{sample}/counts/{sample}.txt.populated.gz",
    log:
        "{folder}/log/plot_mosaic_counts/{sample}.log",
    conda:
        "../../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../../scripts/ashleys/utils/populated_counts_for_qc_plot.py"


rule ashleys_plot_mosaic_counts:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.populated.gz",
        info="{folder}/{sample}/counts/{sample}.info_raw",
    output:
        # "{folder}/{sample}/plots/counts/CountComplete.raw.pdf",
        report(
            "{folder}/{sample}/plots/counts/CountComplete.raw.pdf",
            category="Mosaic Counts",
            subcategory="{sample}",
            labels={"Cell": "ALL", "Type": "raw"},
        ),
    log:
        "{folder}/log/plot_mosaic_counts/{sample}.log",
    params:
        mouse_assembly=True if config["reference"] in ["mm10", "mm39"] else False,
    conda:
        "../../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        """
        LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output}  > {log} 2>&1
        """
