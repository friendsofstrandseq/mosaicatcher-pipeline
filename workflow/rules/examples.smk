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
    shell:
        "unzip {input} -d ."


# TODO: Adapt according reference
rule dl_external_data:
    """
    rule fct: Download External files 
    input: files stored on Zenodo
    output: touch file to check if everything was running correctly
    """
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
            keep_local=True,
        ),
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz",
            keep_local=True,
        ),
    output:
        touch("config_output/dl_external_data.ok"),
    log:
        touch("log/config_output/dl_external_data.ok"),


# TODO: Adapt according reference
rule dl_external_data_index:
    """
    rule fct: Download External files 
    input: files stored on Zenodo
    output: touch file to check if everything was running correctly
    """
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
            keep_local=True,
        ),
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz.tbi",
            keep_local=True,
        ),
    output:
        touch("config_output/dl_external_data_index.ok"),
    log:
        touch("log/config_output/dl_external_data_index.ok"),

#####################


rule download_hg19_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        touch("{output_folder}/config/ref_genomes/download_hg19_reference.ok".format(output_folder=config["output_location"])),
    log:
        touch("{output_folder}/log/ref_genomes/download_hg19_reference.ok".format(output_folder=config["output_location"])),
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg19.fa.gz")

rule download_hg38_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        touch("{output_folder}/config/ref_genomes/download_hg38_reference.ok".format(output_folder=config["output_location"])),
    log:
        touch("{output_folder}/log/ref_genomes/download_hg38_reference.ok".format(output_folder=config["output_location"])),
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg38.fa.gz")

rule download_T2T_reference:
    input:
        HTTP.remote(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        touch("{output_folder}/config/ref_genomes/download_T2T_reference.ok".format(output_folder=config["output_location"])),
    log:
        touch("{output_folder}/log/ref_genomes/download_T2T_reference.ok".format(output_folder=config["output_location"])),
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/T2T.fa.gz")

rule samtools_faindex:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa.gz.fai"
    log:
        "{file}.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools faidx {input}"

