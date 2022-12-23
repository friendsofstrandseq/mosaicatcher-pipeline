# if config["check_sm_tag"] == True:
#     rule check_bam_compliancy:
#         input:
#             bam=lambda wc: expand(
#                 "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
#                 folder=config["data_location"],
#                 sample=wc.sample,
#                 cell=bam_per_sample_local[str(wc.sample)],
#             ),
#         output:
#             "{folder}/{sample}/config/bams_check_sm_rg_tags.tsv"
#         conda:
#             "../envs/mc_base.yaml"
#         script:
#             ""


#     rule change_id_and_sm_tags:
#         input:
#             bam_orig="{path}/{SM}/raw/{ID}.sort.mdup.bam",
#         output:
#             bam_out="{path}/{SM}/bam/{ID}.sort.mdup.bam",
#         conda:
#             "../envs/mc_bioinfo_tools.yaml"
#         resources:
#             mem_mb="16000",
#             time="10:00:00",
#         shell:
#             """
#             # old_id=$(samtools view -H {input.bam_orig} | grep "^@RG" | grep -P -o "\tID:[A-Za-z0-9_\-]*" | sed 's/ID://g')
#             # old_sm=$(samtools view -H {input.bam_orig} | grep "^@RG" | grep -P -o "\tSM:[A-Za-z0-9_\-]*" | sed 's/SM://g')
#             # echo "{wildcards.ID} {wildcards.SM} $old_id $old_sm"
#             # First, the 'ID' tag
#             samtools view -H {input.bam_orig} | sed "s/ID:[A-Za-z0-9_-]*/ID:{wildcards.ID}/g;s/SM:[A-Za-z0-9_-]*/SM:{wildcards.SM}/g;s/RG:Z:[A-Za-z0-9_-]*/RG:Z:{wildcards.ID}/g" > {output.bam_out}.header
#             samtools view -h {input.bam_orig} | sed "s/ID:[A-Za-z0-9_-]*/ID:{wildcards.ID}/g;s/SM:[A-Za-z0-9_-]*/SM:{wildcards.SM}/g;s/RG:Z:[A-Za-z0-9_-]*/RG:Z:{wildcards.ID}/g" | samtools view -bS > {output.bam_out}.core
#             samtools reheader -P {output.bam_out}.header {output.bam_out}.core > {output.bam_out}
#             # Remove intermediate files
#             rm {output.bam_out}.header {output.bam_out}.core
#             """


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
