

rule index_input_bam:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    log:
        "{folder}/log/index_input_bam/{sample}/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule index_haplotag_bam:
    input:
        "{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
    output:
        "{folder}/{sample}/haplotag/bam/{cell}.bam.htg.bai",
    log:
        "{folder}/log/index_haplotag_bam/{sample}/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule compress_indiv_freebayes_vcf:
    input:
        vcf="{folder}/{sample}/snv_genotyping/{chrom}.vcf",
    output:
        vcf="{folder}/{sample}/snv_genotyping/{chrom}.vcf.gz",
    log:
        "{folder}/log/compress_indiv_freebayes_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "bgzip {input.vcf} > {log} 2>&1"


rule index_indiv_freebayes_vcf:
    input:
        vcf="{folder}/{sample}/snv_genotyping/{chrom}.vcf.gz",
    output:
        vcf="{folder}/{sample}/snv_genotyping/{chrom}.vcf.gz.tbi",
    log:
        "{folder}/log/index_indiv_freebayes_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"


rule compress_indiv_strandphaser_vcf:
    input:
        vcf="{folder}/{sample}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf",
    output:
        vcf="{folder}/{sample}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
    log:
        "{folder}/log/compress_indiv_strandphaser_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "bgzip {input.vcf} > {log} 2>&1"


rule index_indiv_strandphaser_vcf:
    input:
        vcf="{folder}/{sample}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz",
    output:
        vcf="{folder}/{sample}/strandphaser/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi",
    log:
        "{folder}/log/index_indiv_strandphaser_vcf/{sample}/{chrom}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"


rule index_merged_strandphaser_vcf:
    input:
        vcf="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz",
    output:
        vcf="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
    log:
        "{folder}/log/index_merged_strandphaser_vcf/phased-snvs/{sample}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "tabix -p vcf {input.vcf} > {log} 2>&1"




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
