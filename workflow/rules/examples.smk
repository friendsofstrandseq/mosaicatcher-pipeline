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
        HTTP.remote("https://sandbox.zenodo.org/record/1074721/files/TEST_EXAMPLE_DATA.zip", keep_local=True)
        # HTTP.remote("https://sandbox.zenodo.org/record/1062186/files/report_TALL.zip", keep_local=True)
    output:
        # directory(config["input_bam_location"])
        touch(config["output_location"] + "config/dl_example_data.ok")
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
        HTTP.remote("https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna", keep_local=True),
        HTTP.remote("https://sandbox.zenodo.org/record/1074721/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz", keep_local=True),
    output:
        touch(config["output_location"] + "config/dl_external_data.ok")

# TODO: Adapt according reference
rule dl_external_data_index:
    """
    rule fct: Download External files 
    input: files stored on Zenodo
    output: touch file to check if everything was running correctly
    """
    input:
        HTTP.remote("https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai", keep_local=True),
        HTTP.remote("https://sandbox.zenodo.org/record/1074721/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz.tbi", keep_local=True),
    output:
        touch(config["output_location"] + "config/dl_external_data_index.ok")
