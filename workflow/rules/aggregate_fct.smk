def aggregate_phased_haps(wildcards):
    df = pd.read_csv(
        checkpoints.summarise_ploidy.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.summary,
        sep="\t",
    )
    df = df.loc[df["50%"] >= 2]
    chrom_list = [e for e in df["#chrom"].values.tolist() if e != "genome"]
    return expand(
        "{{folder}}/{{sample}}/strandphaser/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        chrom=chrom_list,
    )


def aggregate_vcf_gz(wildcards):
    df = pd.read_csv(
        checkpoints.summarise_ploidy.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.summary,
        sep="\t",
    )
    df = df.loc[df["50%"] >= 2]
    chrom_list = [e for e in df["#chrom"].values.tolist() if e != "genome"]
    return expand(
        "{{folder}}/{{sample}}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
        chrom=chrom_list,
    )


def aggregate_vcf_gz_tbi(wildcards):
    df = pd.read_csv(
        checkpoints.summarise_ploidy.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.summary,
        sep="\t",
    )
    df = df.loc[df["50%"] >= 2]
    chrom_list = [e for e in df["#chrom"].values.tolist() if e != "genome"]
    return expand(
        "{{folder}}/{{sample}}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi",
        chrom=chrom_list,
    )


def locate_snv_vcf(wildcards):
    if (
        "snv_calls" not in config["references_data"][config["reference"]]
        or wildcards.sample
        not in config["references_data"][config["reference"]]["snv_calls"]
        or config["references_data"][config["reference"]]["snv_sites_to_genotype"][
            wildcards.sample
        ]
        == ""
    ):
        if (
            "snv_sites_to_genotype" in config["references_data"][config["reference"]]
            and config["references_data"][config["reference"]]["snv_sites_to_genotype"]
            != ""
        ):
            if os.path.isfile(
                config["references_data"][config["reference"]]["snv_sites_to_genotype"]
            ):
                return "{}/{}/snv_genotyping/{}.vcf".format(
                    wildcards.folder, wildcards.sample, wildcards.chrom
                )
            else:
                print("ISSUE")

                return "{}/{}/snv_calls/{}.vcf".format(
                    wildcards.folder, wildcards.sample, wildcards.chrom
                )
        else:
            return "{}/{}/snv_calls/{}.vcf".format(
                wildcards.folder, wildcards.sample, wildcards.chrom
            )
    else:
        return "{}/{}/external_snv_calls/{}.vcf".format(
            wildcards.folder, wildcards.sample, wildcards.chrom
        )


def aggregate_cells_segmentation(wildcards):
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()

    return expand(
        "{folder}/{sample}/segmentation/segmentation-per-cell/{cell}.txt",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )


def aggregate_cells_count_plot(wildcards):
    import pandas as pd

    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )


    # dict_cells_nb_per_sample = df.groupby("sample")["cell"].nunique().to_dict()
    dict_cells_nb_per_sample = {k:len(v) for k,v in cell_per_sample.items()}
    samples = list(dict_cells_nb_per_sample.keys())

    # cell_list = df.cell.tolist()
    cell_list = cell_per_sample[wildcards.sample]
    # tmp_dict = (
    #     df[["sample", "cell"]]
    #     .groupby("sample")["cell"]
    #     .apply(lambda r: sorted(list(r)))
    #     .to_dict()
    # )
    tmp_dict = {
        s: {i + 1: c for i, c in enumerate(cell_list)}
        for s, cell_list in cell_per_sample.items()
    }
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"
    # print(tmp_dict)

    list_indiv_plots = list()

    return [sub_e for e in [ 
        expand(
            "{folder}/{sample}/plots/counts_{plottype}/{cell}.{i}.pdf",
            folder=config["data_location"],
            sample=wildcards.sample,
            plottype=plottype_counts,
            cell=tmp_dict[wildcards.sample][i],
            i=i,
        )
        for i in range(dict_cells_nb_per_sample[wildcards.sample] + 1)
        # for i in range(dict_cells_nb_per_sample[sample] + 1)
    ]for sub_e in e ]





    # cell_list = df.cell.tolist()
    # tmp_dict = (
    #     df[["sample", "cell"]]
    #     .groupby("sample")["cell"]
    #     .apply(lambda r: sorted(list(r)))
    #     .to_dict()
    # )
    # tmp_dict = {
    #     s: {i + 1: c for i, c in enumerate(cell_list)}
    #     for s, cell_list in tmp_dict.items()
    # }
    # for s in tmp_dict.keys():
    #     tmp_dict[s][0] = "SummaryPage"

    # return expand(
    #     "{folder}/{sample}/plots/counts/{cell}.{i}.pdf",
    #     folder=config["data_location"],
    #     sample=wildcards.sample,
    #     cell=cell_list,
    #     i=tmp_dict[i],
    # )
