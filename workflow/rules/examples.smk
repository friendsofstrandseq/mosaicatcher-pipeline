import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule dl_example_data:
    """
    rule fct: Download BAM example data as input for MosaiCatcher pipeline
    input: zip file stored on Zenodo
    output: input_bam_location given by the user
    """
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/TEST_EXAMPLE_DATA.zip",
            keep_local=True,
        ),
    output:
        touch("config_output/dl_example_data.ok"),
    log:
        touch("log/config_output/dl_example_data.ok"),
    conda:
        "../envs/mc_base.yaml"
    run:
        shell("unzip {input} -d .")
        # directory = ".tests/data_example/"
        # if not os.path.exists(directory):
        #     os.makedirs(directory)
        # shell("mv {input} .tests/data_example/TEST_EXAMPLE_DATA.zip")
        # shell("unzip .tests/data_example/TEST_EXAMPLE_DATA.zip")




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
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg19.fa.gz")
        shell("gunzip workflow/data/ref_genomes/hg19.fa.gz")


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
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg38.fa.gz")
        shell("gunzip workflow/data/ref_genomes/hg38.fa.gz")


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
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/T2T.fa.gz")
        shell("gunzip workflow/data/ref_genomes/T2T.fa.gz")


# rule download_T2T_tarball:
#     input:
#         HTTP.remote()
#     output: 
#         ".tar.gz"
#     log:
#     run:

rule install_T2T_tarball:
    input: 
        ".tar.gz"
    output:
        touch("workflow/data/ref_genomes/config/T2T_R_tarball_install.ok")
    log:
        "workflow/data/ref_genomes/log/T2T_R_tarball_install.log"
    conda:
        "../envs/rtools.yaml"
    shell:
        """
        R_path=$(which R | grep -P "\.snakemake" | sed 's/R is //g')
        "$R_path" -e 'install.packages("{input}")  2>&1 > {log}'
        """



rule samtools_faindex:
    input:
        ancient("{file}.fa")
    output:
        "{file}.fa.fai"
    log:
        "{file}.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools faidx {input}"

