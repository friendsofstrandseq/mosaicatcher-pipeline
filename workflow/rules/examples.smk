import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule dl_example_data:
    input:
        # HTTP.remote("https://git.embl.de/tweber/mosaicatcher-update/-/raw/dev/workflow/bam/RPE-BM510/all/BM510x04_PE20301.sort.mdup.bam", keep_local=True)
        HTTP.remote("https://sandbox.zenodo.org/record/1060422/files/input_data.tar.gz", keep_local=True)
    output:
        directory("TEST_EXAMPLE_DATA")
    shell:
        "mkdir  {output};"
        "tar -xf {input} -C {output}"