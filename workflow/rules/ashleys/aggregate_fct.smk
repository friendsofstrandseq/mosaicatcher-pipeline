def aggregate_correct_cells_bam(wildcards):
    if config["use_light_data"] is False:
        df = pd.read_csv(
            checkpoints.tune_predictions_based_on_threshold.get(
                sample=wildcards.sample, folder=config["data_location"]
            ).output[0],
            sep="\t",
        )
        cell_list = df.loc[df["prediction"] == 1].cell.tolist()

        return expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wildcards.sample,
            cell=cell_list,
        )
    else:
        return expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wildcards.sample,
            cell=cell_per_sample[str(wildcards.sample)],
        )


def selected_input_bam(wildcards):
    df = pd.read_csv(
        checkpoints.mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    df = df.loc[df["mapped"] > 0]
    cell_list = df.cell.tolist()
    # # print(cell_list)

    return expand(
        "{folder}/{sample}/bam_ashleys/{cell}.sort.mdup.bam",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )


def aggregate_correct_cells_plot(wildcards):
    if config["use_light_data"] is False:
        df = pd.read_csv(
            checkpoints.tune_predictions_based_on_threshold.get(
                sample=wildcards.sample, folder=config["data_location"]
            ).output[0],
            sep="\t",
        )
        cell_list = df.loc[df["prediction"] == 1].cell.tolist()

        return expand(
            "{folder}/{sample}/plots/alfred/{cell}_gc_{alfred_plot}.png",
            folder=config["data_location"],
            sample=wildcards.sample,
            cell=cell_list,
            alfred_plot=config["alfred_plots"],
        )
    else:
        return expand(
            "{folder}/{sample}/plots/alfred/{cell}_gc_{alfred_plot}.png",
            folder=config["data_location"],
            sample=wildcards.sample,
            cell=cell_per_sample[str(wildcards.sample)],
            alfred_plot=config["alfred_plots"],
        )


def select_binbed(wildcards):
    return get_bin_bed_file()


def select_ashleys_labels(wildcards):
    # if bypass_ashleys is True > pick labels_ashleys_bypass.tsv
    if config["bypass_ashleys"] or config["use_light_data"]:
        return expand(
            "{folder}/{sample}/cell_selection/labels_ashleys_bypass.tsv",
            folder=config["data_location"],
            sample=wildcards.sample,
        )
    else:
        # if bypass_ashleys is False > pick labels_ashleys.tsv
        return expand(
            "{folder}/{sample}/cell_selection/labels_ashleys.tsv",
            folder=config["data_location"],
            sample=wildcards.sample,
        )
