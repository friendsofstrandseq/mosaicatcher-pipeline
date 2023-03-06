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


rule install_T2T_BSgenome_tarball:
    input:
        tarball="workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
    output:
        touch("workflow/data/ref_genomes/config/T2T_R_tarball_install.ok"),
    log:
        "workflow/data/ref_genomes/log/T2T_R_tarball_install.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/utils/install_R_tarball.R"


rule empty_install:
    output:
        touch("workflow/data/ref_genomes/config/fake_install.ok"),
    log:
        "workflow/data/ref_genomes/config/fake_install.log",


rule samtools_faindex:
    input:
        ancient("{file}.fa"),
    output:
        "{file}.fa.fai",
    log:
        "{file}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        "samtools faidx {input}"


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
        touch("workflow/data/scNOVA/log/dl.ok"),
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/scNOVA_data_models.zip
        unzip workflow/data/scNOVA_data_models.zip
        """
