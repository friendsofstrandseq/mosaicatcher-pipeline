# from workflow.scripts.utils.utils import get_mem_mb

################################################################################
# UTILS                                                                        #
################################################################################


<<<<<<< HEAD
rule index_input_bam:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{input_folder}/{sample}/all/{cell}.sort.mdup.bam",
    output:
        "{input_folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
    log:
        "{input_folder}/log/index_input_bam/{sample}/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    # resources:
    #     mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule index_haplotag_bam:
=======


rule index_bam:
>>>>>>> master
    """
    rule fct:
    input:
    output:
    """
    input:
<<<<<<< HEAD
        "{output_folder}/haplotag/bam/{sample}/{cell}.bam",
=======
        ancient("{file}.bam")
>>>>>>> master
    output:
        "{output_folder}/haplotag/bam/{sample}/{cell}.bam.bai",
    log:
        "{output_folder}/log/index_haplotag_bam/{sample}/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    # resources:
    #     mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule index_merged_bam:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{output_folder}/merged_bam/{sample}/merged.bam",
    output:
        "{output_folder}/merged_bam/{sample}/merged.bam.bai",
    log:
        "{output_folder}/log/merged_bam/{sample}/merged.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    # resources:
    #     mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule compress_indiv_freebayes_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{output_folder}/snv_genotyping/{sample}/{chrom}.vcf",
    output:
        vcf="{output_folder}/snv_genotyping/{sample}/{chrom}.vcf.gz",
    log:
        "{output_folder}/log/compress_indiv_freebayes_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "bgzip {input.vcf} > {log} 2>&1"


rule index_indiv_freebayes_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
<<<<<<< HEAD
        vcf="{output_folder}/snv_genotyping/{sample}/{chrom}.vcf.gz",
=======
        vcf=ancient("{file}.vcf"),
>>>>>>> master
    output:
        vcf="{output_folder}/snv_genotyping/{sample}/{chrom}.vcf.gz.tbi",
    log:
        "{output_folder}/log/index_indiv_freebayes_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"


rule compress_indiv_strandphaser_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{output_folder}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf",
    output:
        vcf="{output_folder}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
    log:
        "{output_folder}/log/compress_indiv_strandphaser_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "bgzip {input.vcf} > {log} 2>&1"


rule index_indiv_strandphaser_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
<<<<<<< HEAD
        vcf="{output_folder}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
=======
        vcf=ancient("{file}.vcf.gz"),
>>>>>>> master
    output:
        vcf="{output_folder}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi",
    log:
        "{output_folder}/log/index_indiv_strandphaser_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"


rule index_merged_strandphaser_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{output_folder}/strandphaser/phased-snvs/{sample}.vcf.gz",
    output:
        vcf="{output_folder}/strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
    log:
        "{output_folder}/log/index_merged_strandphaser_vcf/phased-snvs/{sample}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"
