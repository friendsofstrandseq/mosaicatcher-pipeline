import math


rule segmentation:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
    output:
        "{folder}/{sample}/segmentation/{sample}.txt.fixme",
    log:
        "{folder}/log/segmentation/{sample}/{sample}.log",
    params:
        min_num_segs=lambda wc: math.ceil(200000 / float(config["window"])),  # bins to represent 200 kb
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input.counts} > {log} 2>&1
        """


# FIXME: This is a workaround because latest versions of "mosaic segment" don't compute the "bps" column properly. Remove once fixed in the C++ code.
rule fix_segmentation:
    input:
        ancient("{folder}/{sample}/segmentation/{sample}.txt.fixme"),
    output:
        "{folder}/{sample}/segmentation/{sample}.txt",
    log:
        "{folder}/log/segmentation/{sample}/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    params:
        script="workflow/scripts/segmentation_scripts/fix_segmentation.awk",
        window=config["window"],
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v name={wildcards.sample} -v window={params.window} -f {params.script} {input} > {output}
        """


rule segment_one_cell:
    input:
        "{folder}/{sample}/counts/counts-per-cell/{cell}.txt.percell.gz",
    output:
        "{folder}/{sample}/segmentation/segmentation-per-cell/{cell}.txt",
    log:
        "{folder}/log/segmentation/{sample}/segmentation-per-cell/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    params:
        min_num_segs=lambda wc: math.ceil(200000 / float(config["window"])),  # bins to represent 200 kb
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher segment \
        --remove-none \
        --forbid-small-segments {params.min_num_segs} \
        -M 50000000 \
        -o {output} \
        {input} > {log} 2>&1
        """


rule segmentation_selection:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.gz",
        jointseg="{folder}/{sample}/segmentation/{sample}.txt",
        singleseg=aggregate_cells_segmentation,
        info="{folder}/{sample}/counts/{sample}.info",
    output:
        jointseg="{folder}/{sample}/segmentation/Selection_jointseg.txt",
        singleseg="{folder}/{sample}/segmentation/Selection_singleseg.txt",
        strand_states="{folder}/{sample}/segmentation/Selection_initial_strand_state",
    log:
        "{folder}/log/segmentation/segmentation_selection/{sample}.log",
    params:
        cellnames=lambda wc: ",".join(
            [
                cell
                for cell in cell_per_sample[wc.sample]
                if cell
                in [
                    e.split(config["abs_path"])[-1].split(".")[0]
                    for e in aggregate_cells_segmentation(wc)
                ]
            ]
        ),
        sce_min_distance=config["sce_min_distance"],
        additional_sce_cutoff=config["additional_sce_cutoff"],
        min_diff_jointseg=config["min_diff_jointseg"],
        min_diff_singleseg=config["min_diff_singleseg"],
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        PYTHONPATH="" # Issue #1031 (https://bitbucket.org/snakemake/snakemake/issues/1031)
        python workflow/scripts/segmentation_scripts/detect_strand_states.py \
            --sce_min_distance {params.sce_min_distance} \
            --sce_add_cutoff {params.additional_sce_cutoff} \
            --min_diff_jointseg {params.min_diff_jointseg} \
            --min_diff_singleseg {params.min_diff_singleseg} \
            --output_jointseg {output.jointseg} \
            --output_singleseg {output.singleseg} \
            --output_strand_states {output.strand_states} \
            --samplename {wildcards.sample} \
            --cellnames {params.cellnames} \
            {input.info} \
            {input.counts} \
            {input.jointseg} \
            {input.singleseg} > {log} 2>&1
        """
