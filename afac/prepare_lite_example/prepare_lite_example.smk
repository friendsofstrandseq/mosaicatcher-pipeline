REF = "workflow/data/ref_genomes/hg38.fa"
folder = "/g/korbel2/weber/workspace/mosaicatcher-update/TEST_EXAMPLE_DATA/RPE-BM510"
chrom = ["chr17"]

import os

list_dir = [
    e.split("/")[-1].replace(".sort.mdup.bam", "")
    for e in os.listdir("{folder}/all/".format(folder=folder))
    if e.endswith(".sort.mdup.bam")
]
print(list_dir)
print(
    expand(
        "{folder}/lite_final_{chrom}/{cell}.sort.mdup.bam",
        folder=folder,
        cell=list_dir,
        chrom=chrom,
    )
)


rule all:
    input:
        expand(
            "{folder}/lite_final_{chrom}/{cell}.sort.mdup.bam.bai",
            folder=folder,
            cell=list_dir,
            chrom=chrom,
        ),


rule extract_chr_reads:
    input:
        "{folder}/{sample}/all/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/filter_{chr}/{cell}.bam",
    shell:
        "samtools view {input} {wildcards.chr} -b | samtools sort -n -o {output}"


rule bgzip:
    input:
        "{file}.fastq",
    output:
        "{file}.fastq.gz",
    shell:
        "bgzip {input}"


rule bam_to_fastq:
    input:
        "{folder}/{sample}/filter_{chr}/{cell}.bam",
    output:
        fq1="{folder}/{sample}/fastq_{chr}/{cell}.{chr}.1.fastq",
        fq2="{folder}/{sample}/fastq_{chr}/{cell}.{chr}.2.fastq",
    log:
        "{folder}/{sample}/log/fastq_{chr}/{cell}.{chr}.2.log",
    shell:
        "bamToFastq -i {input} -fq {output.fq1} -fq2 {output.fq2} 2>&1 > {log}"
        # "bamToFastq -i {input} -fq {output.fq1} 2>&1 > {log}"


rule remap:
    input:
        fq1="{folder}/{sample}/fastq_{chr}/{cell}.{chr}.1.fastq.gz",
        fq2="{folder}/{sample}/fastq_{chr}/{cell}.{chr}.2.fastq.gz",
    output:
        "{folder}/{sample}/remap_{chr}/{cell}.sort.bam",
        # "{folder}/{sample}/remap_{chr}/{cell}.sor t.mdup.bam"
    params:
        # sample = lambda wc: str(wc.folder).split('/')[-1],
        sample="RPE-BM510",
        REF=REF,
    shell:
        'bwa mem -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:RPE-BM510" {params.REF} {input.fq1} {input.fq2} | samtools sort | samtools view -b -o {output}'
        # 'bwa mem -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:RPE-BM510" {params.REF} {input.fq1} | samtools sort | samtools view -b -o {output}'
        # 'bwa mem -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{params.sample}" {params.REF} {input.fq1} {input.fq2} | samtools sort | samtools view -b -o {output}'


rule markdup:
    input:
        "{folder}/{sample}/remap_{chr}/{cell}.sort.bam",
    output:
        "{folder}/{sample}/lite_final_{chr}/{cell}.sort.mdup.bam",
    conda:
        "sambamba"
    shell:
        "sambamba markdup {input} {output}"


rule index_bam:
    input:
        "{folder}/{sample}/lite_final_{chr}/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/lite_final_{chr}/{cell}.sort.mdup.bam.bai",
    shell:
        "samtools index {input}"
