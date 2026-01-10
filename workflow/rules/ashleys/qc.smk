rule ashleys_generate_features:
    input:
        bam=selected_input_bam,
    output:
        "{folder}/{sample}/predictions/ashleys_features.tsv",
    log:
        "{folder}/log/ashleys/{sample}/features.log",
    conda:
        "../../envs/mc_base.yaml"
    threads: 64
    params:
        windows="5000000 2000000 1000000 800000 600000 400000 200000",
        extension=".sort.mdup.bam",
        folder=lambda wildcards, input: "{}bam_ashleys".format(
            input.bam[0].split("bam_ashleys")[0]
        ),
    resources:
        mem_mb=get_mem_mb_heavy,
        runtime=3600,
    shell:
        "ashleys -j {threads} features -f {params.folder} -w {params.windows} -o {output} --recursive_collect -e {params.extension}"


rule ashleys_predict:
    input:
        folder="{folder}/{sample}/predictions/ashleys_features.tsv",
    output:
        "{folder}/{sample}/cell_selection/labels_ashleys.tsv",
    log:
        "{folder}/log/ashleys/{sample}/prediction_ashleys.log",
    conda:
        "../../envs/mc_base.yaml"
    params:
        model_default="./workflow/ashleys_models/svc_default.pkl",
        model_stringent="./workflow/ashleys_models/svc_stringent.pkl",
    resources:
        mem_mb=get_mem_mb,
        runtime=3600,
    shell:
        "ashleys predict -p {input.folder} -o {output} -m {params.model_default}"


rule ashleys_generate_default_labels:
    input:
        bam=selected_input_bam,
    output:
        "{folder}/{sample}/cell_selection/labels_ashleys_bypass.tsv",
    log:
        "{folder}/log/generate_default_labels/{sample}.log",
    conda:
        "../../envs/mc_base.yaml"
    shell:
        """
        echo "cell\tprediction\tprobability\tsample" > {output}
        for bam in {input.bam} ; do
            # remove path only
            cell=$(basename $bam)
            echo -e "$cell\t1\t1\t{wildcards.sample}" >> {output}
        done
        """


if config["use_light_data"] is False:

    rule ashleys_positive_negative_control_bypass:
        input:
            labels=select_ashleys_labels,
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            labels_corrected="{folder}/{sample}/cell_selection/labels_positive_control_corrected.tsv",
            bypass_cell="{folder}/{sample}/config/bypass_cell.txt",
        log:
            "{folder}/log/positive_control_bypass/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        script:
            "../../scripts/ashleys/utils/positive_negative_control_bypass.py"

    checkpoint tune_predictions_based_on_threshold:
        input:
            "{folder}/{sample}/cell_selection/labels_positive_control_corrected.tsv",
        output:
            "{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/cp_predictions/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        script:
            "../../scripts/ashleys/utils/tune_predictions_based_on_threshold.py"

elif config["use_light_data"] is True:

    rule ashleys_dev_all_cells_correct:
        input:
            # folder="{folder}/{sample}/cell_selection/labels_notebook.tsv",
            folder=select_ashleys_labels,
        output:
            folder="{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/dev_all_cells_correct/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        script:
            "../../scripts/ashleys/utils/dev_all_cells_correct.py"


rule ashleys_plot_plate:
    input:
        labels="{folder}/{sample}/cell_selection/labels.tsv",
    output:
        predictions=report(
            "{folder}/{sample}/plots/plate/ashleys_plate_predictions.pdf",
            category="Ashleys plate plots",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Plot Type": "Predictions"},
        ),
        probabilities=report(
            "{folder}/{sample}/plots/plate/ashleys_plate_probabilities.pdf",
            category="Ashleys plate plots",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Plot Type": "Probabilities"},
        ),
        well_table="{folder}/{sample}/plots/plate/ashleys_well_table.tsv",
    log:
        "{folder}/log/plot_plate/{sample}.log",
    conda:
        "../../envs/rtools.yaml"
    script:
        "../../scripts/ashleys/plotting/plot_plate.R"


if config["publishdir"] != "":

    rule ashleys_publishdir_outputs_ashleys:
        input:
            list_publishdir=publishdir_fct,
        output:
            touch("{folder}/{sample}/config/publishdir_outputs_ashleys.ok"),
        log:
            "{folder}/log/publishdir_outputs_ashleys/{sample}.log",
        conda:
            "../../envs/mc_base.yaml"
        script:
            "../../scripts/ashleys/utils/publishdir.py"


rule ashleys_save_config:
    input:
        "config/config.yaml",
    output:
        "{folder}/{sample}/config/config_ashleys.yaml",
    log:
        "{folder}/log/save_config/{sample}.log",
    conda:
        "../../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../../scripts/ashleys/utils/dump_config.py"


def get_ashleys_outputs(wildcards):
    """
    Get final outputs for ashleys pipeline
    """
    final_list = []

    # Cell selection labels
    final_list.append(
        "{folder}/{sample}/cell_selection/labels.tsv".format(
            folder=config["data_location"], sample=wildcards.sample
        )
    )

    # Config
    final_list.append(
        "{folder}/{sample}/config/config_ashleys.yaml".format(
            folder=config["data_location"], sample=wildcards.sample
        )
    )

    # MultiQC if enabled
    if config["MultiQC"] is True:
        final_list.append(
            "{folder}/{sample}/multiqc/multiqc_report/multiqc_report.html".format(
                folder=config["data_location"], sample=wildcards.sample
            )
        )

    # Plate plots if applicable
    sample = wildcards.sample
    if len(cell_per_sample[sample]) in [96, 384]:
        if config["use_light_data"] is False:
            final_list.extend(
                [
                    "{folder}/{sample}/plots/plate/ashleys_plate_predictions.pdf".format(
                        folder=config["data_location"], sample=sample
                    ),
                    "{folder}/{sample}/plots/plate/ashleys_plate_probabilities.pdf".format(
                        folder=config["data_location"], sample=sample
                    ),
                ]
            )

    # Publishdir if configured
    if config["publishdir"] != "":
        final_list.append(
            "{folder}/{sample}/config/publishdir_outputs_ashleys.ok".format(
                folder=config["data_location"], sample=wildcards.sample
            )
        )

    return final_list


rule ashleys_all:
    input:
        get_ashleys_outputs,
    default_target: True
