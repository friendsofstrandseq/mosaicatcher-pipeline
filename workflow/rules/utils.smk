# from workflow.scripts.utils.utils import get_mem_mb

################################################################################
# UTILS                                                                        #
################################################################################


rule index_bam:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    log:
        "{file}.bam.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule compress_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{file}.vcf",
    output:
        vcf="{file}.vcf.gz",
    log:
        "{file}.vcf.gz.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "bgzip {input.vcf} > {log} 2>&1"


rule index_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{file}.vcf.gz",
    output:
        vcf="{file}.vcf.gz.tbi",
    log:
        "{file}.vcf.gz.tbi.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"
