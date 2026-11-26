# Cell selection with ashleys-qc and/or EVAs (docs/cell_selection.md).
#
# The settings come from the cell_selection config block, parsed and validated
# in common.smk (CELL_SELECTION). Where labels.tsv is written:
#   - FASTQ input (ashleys_pipeline: True): checkpoint
#     tune_predictions_based_on_threshold (rules/ashleys/qc.smk)
#   - BAM input with tools: [evas]: rule cell_selection_decision below
#   - BAM input otherwise: touch_labels (count.smk), every cell kept, as before
# Both decision steps run scripts/cell_selection/decide_cell_selection.py and
# also write labels_detailed.tsv (per-tool columns), used by the plate plot.


if "evas" in CELL_SELECTION["tools"]:

    # EVAs is not on bioconda yet: it runs from a conda environment installed by
    # hand (cell_selection.evas_env), or from the PATH when evas_env is empty.
    # It runs outside the pipeline container, which does not include it.
    # EVAs reads {parent-dir}/{sample}/counts/{sample}.txt.raw.gz, i.e. the raw
    # mosaicatcher counts, and a 200 kb GC matrix with a "GC" column: the
    # pipeline GC matrix is passed with its "%GC" column renamed.
    rule evas_predict:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
            gc=ancient(get_gc_matrix_file()),
        output:
            predictions="{folder}/{sample}/cell_selection/evas/EVAs_{sample}/{sample}_predictions_with_reasons_and_bins.tsv",
        log:
            "{folder}/log/evas_predict/{sample}.log",
        params:
            out_dir=lambda wc: f"{wc.folder}/{wc.sample}/cell_selection/evas",
            extra=CELL_SELECTION["evas_extra_args"],
        container:
            None
        conda:
            CELL_SELECTION["evas_env"] or None
        resources:
            mem_mb=get_mem_mb_heavy,
            runtime=60,
        shell:
            """
            exec > {log} 2>&1
            if ! command -v evas > /dev/null 2>&1; then
                echo "ERROR: 'evas' command not found. EVAs is not on bioconda yet:"
                echo "install it by hand (docs/cell_selection.md) and set"
                echo "cell_selection: evas_env: /path/to/conda/env"
                exit 1
            fi
            mkdir -p {params.out_dir}
            gc={params.out_dir}/gc_matrix_200000.tsv
            zcat -f {input.gc} | sed '1s/%GC/GC/' > $gc
            evas \
                --sample {wildcards.sample} \
                --parent-dir {wildcards.folder} \
                --out-dir {params.out_dir} \
                --gc-fname $gc \
                {params.extra}
            rm -f $gc
            """


if config["ashleys_pipeline"] is False and "evas" in CELL_SELECTION["tools"]:

    rule cell_selection_decision:
        input:
            unpack(cell_selection_inputs),
        output:
            labels="{folder}/{sample}/cell_selection/labels.tsv",
            detailed="{folder}/{sample}/cell_selection/labels_detailed.tsv",
        log:
            "{folder}/log/cell_selection_decision/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            decision=CELL_SELECTION["decision"],
            ashleys_threshold=config["ashleys_threshold"],
            evas_borderline_as_pass=CELL_SELECTION["evas_borderline_as_pass"],
        script:
            "../scripts/cell_selection/decide_cell_selection.py"


if cell_selection_active():

    rule plot_plate_cell_selection:
        input:
            detailed="{folder}/{sample}/cell_selection/labels_detailed.tsv",
        output:
            pdf=report(
                "{folder}/{sample}/plots/plate/cell_selection_plate.pdf",
                category="Cell selection",
                subcategory="{sample}",
                labels={"Sample": "{sample}", "Plot Type": "Plate"},
            ),
            well_table="{folder}/{sample}/plots/plate/cell_selection_well_table.tsv",
        log:
            "{folder}/log/plot_plate_cell_selection/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            ashleys_threshold=config["ashleys_threshold"],
        script:
            "../scripts/cell_selection/plot_plate_cell_selection.py"
