def aggregate_phased_haps(wildcards):
    """
    Function based on checkpoint summarise_ploidy to process only chromosomes where
    the median ploidy status is equal or above 2 for all segments
    Return phased_haps.txt as input for combine_strandphaser_output
    """
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
    """
    Function based on checkpoint summarise_ploidy to process only chromosomes where
    the median ploidy status is equal or above 2 for all segments
    Return {chrom}_phased.vcf.gz as input for merge_strandphaser_vcfs
    """
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
    """
    Function based on checkpoint summarise_ploidy to process only chromosomes where
    the median ploidy status is equal or above 2 for all segments
    Return {chrom}_phased.vcf.gz.tbi as input for merge_strandphaser_vcfs
    """
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
    """
    Function as an input in run_strandphaser_per_chrom
    Trigger based on config file either:
        - regenotyping of 1000G file / other vcf file using freebayes
        - de novo calling using bcftools
    """
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
                # print("ISSUE")

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
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
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


def aggregate_cells_haplotag_tables(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()
    # print(cell_list)
    # print(expand(
    #     "{folder}/{sample}/haplotag/table/by-cell/{cell}.tsv",
    #     folder=config["data_location"],
    #     sample=wildcards.sample,
    #     cell=cell_list,
    # ))
    return expand(
        "{folder}/{sample}/haplotag/table/by-cell/{cell}.tsv",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )


def aggregate_cells_scTRIP_multiplot(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()

    return expand(
        "{folder}/{sample}/plots/scTRIP_multiplot/{cell}/{chrom}.png",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
        chrom=config["chromosomes"],
    )


def unselected_input_bam(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info_removed,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()
    # # print(cell_list)

    # if len(cell_list)>0:
    return expand(
        "{folder}/{sample}/selected/{cell}.sort.mdup.bam",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )
    # else:
    #     return ""


def unselected_input_bai(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info_removed,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()

    # if len(cell_list)>0:
    return expand(
        "{folder}/{sample}/selected/{cell}.sort.mdup.bam.bai",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )
    # else:
    #     return ""


def selected_input_bam(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()
    # # print(cell_list)

    return expand(
        "{folder}/{sample}/selected/{cell}.sort.mdup.bam",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )


def selected_input_bai(wildcards):
    """
    Function based on checkpoint filter_bad_cells_from_mosaic_count
    to process the segmentation only on cells that were flagged as high-quality
    Return {cell}.txt
    """
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()

    return expand(
        "{folder}/{sample}/selected/{cell}.sort.mdup.bam.bai",
        folder=config["data_location"],
        sample=wildcards.sample,
        cell=cell_list,
    )


def remove_unselected_fct(wildcards):
    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, folder=config["data_location"]
        ).output.info_removed,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()
    # print(cell_list)
    # print(len(cell_list))
    if len(cell_list) == 0:
        # if len(cell_list) == 0 or config["input_bam_legacy"] is True:
        return "{folder}/{sample}/config/remove_unselected_bam_empty.ok"
    else:
        return "{folder}/{sample}/config/remove_unselected_bam.ok"


def select_counts_for_SV_calling(wildcards):
    if config["multistep_normalisation_for_SV_calling"] == True:
        return "{folder}/{sample}/counts/multistep_normalisation/{sample}.txt.scaled.GC.VST.reformat.gz"
    else:
        return "{folder}/{sample}/counts/{sample}.txt.raw.gz"


def bsgenome_install(wildcards):
    if config["reference"] == "T2T":
        return "workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz"
    else:
        return "workflow/data/ref_genomes/log/fake_package.ok"


def select_binbed(wildcards):
    if config["reference"] != "mm10":
        return "workflow/data/bin_200kb_all.bed"
    else:
        return "workflow/data/mm10.bin_200kb_all.bed"


def select_labels(wildcards):
    if config["use_strandscape_labels"]:
        return "{folder}/{sample}/cell_selection/labels_strandscape.tsv"
    else:
        return "{folder}/{sample}/cell_selection/labels.tsv"
