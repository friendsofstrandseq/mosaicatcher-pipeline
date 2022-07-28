def aggregate_phased_haps(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output_folder=config["output_location"]
    ).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split("\t")[1]
        if sex == "M":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrX", "chrY"]
            ]
        elif sex == "F":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrY"]
            ]
        return expand(
            "{{output_folder}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
            chrom=config["chromosomes"],
        )


def aggregate_vcf_gz(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output_folder=config["output_location"]
    ).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split("\t")[1]
        if sex == "M":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrX", "chrY"]
            ]
        elif sex == "F":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrY"]
            ]
        return expand(
            "{{output_folder}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
            chrom=config["chromosomes"],
        )


def aggregate_vcf_gz_tbi(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output_folder=config["output_location"]
    ).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split("\t")[1]
        if sex == "M":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrX", "chrY"]
            ]
        elif sex == "F":
            config["chromosomes"] = [
                c for c in config["chromosomes"] if c not in ["chrY"]
            ]
        return expand(
            "{{output_folder}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi",
            chrom=config["chromosomes"],
        )


def aggregate_cells_segmentation(wildcards):
    import pandas as pd

    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, output_folder=config["output_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )
    cell_list = df.cell.tolist()

    return [
        sub_e
        for e in [
            expand(
                "{output_folder}/segmentation/{sample}/segmentation-per-cell/{cell}.txt",
                output_folder=config["output_location"],
                sample=samples,
                cell=cell_list,
            )
            for sample in samples
        ]
        for sub_e in e
    ]


def aggregate_cells_count_plot(wildcards):
    import pandas as pd

    df = pd.read_csv(
        checkpoints.filter_bad_cells_from_mosaic_count.get(
            sample=wildcards.sample, output_folder=config["output_location"]
        ).output.info,
        skiprows=13,
        sep="\t",
    )

    cell_list = df.cell.tolist()
    tmp_dict = (
        df[["sample", "cell"]]
        .groupby("sample")["cell"]
        .apply(lambda r: sorted(list(r)))
        .to_dict()
    )
    tmp_dict = {
        s: {i + 1: c for i, c in enumerate(cell_list)}
        for s, cell_list in tmp_dict.items()
    }
    for s in tmp_dict.keys():
        tmp_dict[s][0] = "SummaryPage"

    return [
        sub_e
        for e in [
            expand(
                "{output_folder}/plots/{sample}/counts/{cell}.{i}.pdf",
                output_folder=config["output_location"],
                sample=samples,
                cell=cell_list,
                i=tmp_dict[i],
            )
            for sample in samples
        ]
        for sub_e in e
    ]
