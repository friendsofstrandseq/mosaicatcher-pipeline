import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule dl_example_data:
    input:
        # HTTP.remote("https://git.embl.de/tweber/mosaicatcher-update/-/raw/dev/workflow/bam/RPE-BM510/all/BM510x04_PE20301.sort.mdup.bam", keep_local=True)
<<<<<<< HEAD
        HTTP.remote("https://sandbox.zenodo.org/record/1060987/files/TEST_EXAMPLE_DATA.zip", keep_local=True)
=======
        HTTP.remote("https://sandbox.zenodo.org/record/1060422/files/input_data.tar.gz", keep_local=True)
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
    output:
        directory(config["input_bam_location"])
    shell:
        "mkdir  {output};"
<<<<<<< HEAD
        "tar -xf {input} -C .;"
        "mv TEST_EXAMPLE_DATA/ {output}"

rule dl_external_data:
    input:
        # HTTP.remote("https://git.embl.de/tweber/mosaicatcher-update/-/raw/dev/workflow/bam/RPE-BM510/all/BM510x04_PE20301.sort.mdup.bam", keep_local=True)
        HTTP.remote("https://sandbox.zenodo.org/record/1060987/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna", keep_local=True),
        HTTP.remote("https://sandbox.zenodo.org/record/1060987/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai", keep_local=True),
        HTTP.remote("https://sandbox.zenodo.org/record/1060987/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz", keep_local=True),
        HTTP.remote("https://sandbox.zenodo.org/record/1060987/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz.tbi", keep_local=True),
    output:
        "sandbox.zenodo.org/record/1060987/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        "sandbox.zenodo.org/record/1060987/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
        "sandbox.zenodo.org/record/1060987/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz"
        "sandbox.zenodo.org/record/1060987/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz.tbi"
    shell:
        "echo Download completed"
=======
        "tar -xf {input} -C {output}"
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
