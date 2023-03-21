import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule dl_example_data:
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/TEST_EXAMPLE_DATA.zip",
            keep_local=True,
        ),
    output:
        touch("config/dl_example_data.ok"),
    log:
        touch("log/config/dl_example_data.ok"),
    run:
        shell("unzip {input} -d .")


rule download_hg19_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg19.fa",
    log:
        "workflow/data/ref_genomes/log/hg19.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg19.fa.gz
        gunzip workflow/data/ref_genomes/hg19.fa.gz
        """


rule download_hg38_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg38.fa",
    log:
        "workflow/data/ref_genomes/log/hg38.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg38.fa.gz
        gunzip workflow/data/ref_genomes/hg38.fa.gz
        """


rule download_T2T_reference:
    input:
        HTTP.remote(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/T2T.fa",
    log:
        "workflow/data/ref_genomes/log/T2T.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/T2T.fa.gz
        gunzip workflow/data/ref_genomes/T2T.fa.gz
        """


rule download_T2T_tarball:
    input:
        HTTP.remote(
            "https://zenodo.org/record/7697400/files/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
    log:
        "workflow/data/ref_genomes/log/T2T_tarball.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz
        """



rule download_arbigent_mappability_track:
    input:
        HTTP.remote(
            "https://zenodo.org/record/7697400/files/mapping_counts_allchrs_hg38.txt",
            keep_local=True,
        ),
    output:
        config["arbigent_data"]["arbigent_mapability_track"],
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    shell:
        """
        directory="workflow/data/arbigent/"
        mkdir -p "$directory"
        mv {input} {output}
        """


rule download_scnova_data:
    input:
        HTTP.remote(
            "https://zenodo.org/record/7697400/files/scNOVA_data_models.zip",
            keep_local=True,
        ),
    output:
        "workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    conda:
        "../envs/scNOVA/scNOVA_DL.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/scNOVA_data_models.zip
        unzip workflow/data/scNOVA_data_models.zip -d workflow/data/
        """
