def aggregate_phased_haps(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output=config["output_location"]
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
            "{{output}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
            chrom=config["chromosomes"],
        )


def aggregate_vcf_gz(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output=config["output_location"]
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
            "{{output}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
            chrom=config["chromosomes"],
        )


def aggregate_vcf_gz_tbi(wildcards):
    with checkpoints.determine_sex_per_cell.get(
        sample=wildcards.sample, output=config["output_location"]
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
            "{{output}}/strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi",
            chrom=config["chromosomes"],
        )
