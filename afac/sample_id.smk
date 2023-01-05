folder = "/g/korbel/weber/MosaiCatcher_files/POOLING/PSEUDOPOOL/pseudopool"

import os

list_dir = [
    e.split("/")[-1].replace(".sort.mdup.bam", "")
    for e in os.listdir("{folder}/all/".format(folder=folder))
    if e.endswith(".sort.mdup.bam")
]


rule all:
    input:
        expand(
            "{folder}/vcf/{file}.vcf",
            folder=folder,
            file=list_dir[1],
        ),


rule freebayes:
    input:
        samples="{folder}/all/{file}.sort.mdup.bam",
        sites="/g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz",
        ref="/g/korbel2/weber/workspace/mosaicatcher-update/workflow/data/ref_genomes/hg38.fa",
    output:
        vcf="{folder}/vcf/{file}.vcf",
    params:
        extra="-@ {input.sites} --only-use-input-alleles --genotype-qualities",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 12
    wrapper:
        "v1.12.2/bio/freebayes"


# rule regenotype_SNVs:
#     input:
#         bam="{folder}/all/{file}.sort.mdup.bam",
#         sites="/g/korbel2/weber/MosaiCatcher_files/EXTERNAL_DATA/snv_sites_to_genotype/ALL.chr1-22plusXY_T2T_sites.20170504.renamedCHR.vcf.gz",
#         fasta="/g/korbel2/weber/workspace/mosaicatcher-update/workflow/data/ref_genomes/hg38.fa",
#     output:
#         vcf="{folder}/vcf/{file}.vcf",
#     resources:
#         mem_mb="8G",
#         time="10:00:00",
#     conda:
#         "../workflow/envs/mc_bioinfo_tools.yaml"
#     shell:
#         """
#         freebayes \
#             -f {input.fasta} \
#             -@ {input.sites} \
#             --only-use-input-alleles {input.bam} \
#             --genotype-qualities \
#         > {output.vcf}
#         """
# """
# (freebayes \
#     -f {input.fasta} \
#     -r {wildcards.chrom} \
#     -@ {input.sites} \
#     --only-use-input-alleles {input.bam} \
#     --genotype-qualities \
# | bcftools view \
#     --exclude-uncalled \
#     --types snps \
#     --genotype het \
#     --include "QUAL>=10" \
# > {output.vcf}) 2> {log}
# """
