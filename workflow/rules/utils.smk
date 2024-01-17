rule check_sm_tag:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/checks/{cell}.sm_check.ok",
    log:
        "{folder}/log/{sample}/checks/{cell}.sm_check.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        """
        sample_name="{wildcards.sample}"
        sm_tag=$(samtools view -H {input} | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g")

        if [[ $sample_name == $sm_tag ]]; then
            echo "{input}: $sm_tag $sample_name OK" > {output}
            echo "{input}: $sm_tag $sample_name OK" > {log}
        else
            echo "{input}: $sm_tag $sample_name MISMATCH" > {log}
        fi
        """


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


########################################


rule save_config:
    input:
        "config/config.yaml",
    output:
        "{folder}/{sample}/config/config.yaml",
    log:
        "{folder}/log/save_config/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/utils/dump_config.py"


rule save_conda_versions_mc_base:
    output:
        "{folder}/{sample}/config/conda_export/mc_base.yaml",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "conda env export > {output}"


rule save_conda_versions_mc_bioinfo_tools:
    output:
        "{folder}/{sample}/config/conda_export/mc_bioinfo_tools.yaml",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "conda env export > {output}"


rule save_conda_versions_rtools:
    output:
        "{folder}/{sample}/config/conda_export/rtools.yaml",
    conda:
        "../envs/rtools.yaml"
    shell:
        "conda env export > {output}"


if config["scNOVA"] == True:

    rule save_conda_versions_scNOVA_bioinfo_tools:
        output:
            "{folder}/{sample}/config/conda_export/scNOVA.yaml",
        conda:
            "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
        shell:
            "conda env export > {output}"

    rule save_conda_versions_scNOVA_DL:
        output:
            "{folder}/{sample}/config/conda_export/scNOVA_DL.yaml",
        conda:
            "../envs/scNOVA/scNOVA_DL.yaml"
        shell:
            "conda env export > {output}"

    rule save_conda_versions_scNOVA_R:
        output:
            "{folder}/{sample}/config/conda_export/scNOVA_R.yaml",
        conda:
            "../envs/scNOVA/scNOVA_R.yaml"
        shell:
            "conda env export > {output}"


# PUBLISHDIR

if config["publishdir"] != "":

    rule publishdir_outputs_mc:
        input:
            list_publishdir=publishdir_fct_mc,
        output:
            touch("{folder}/{sample}/config/publishdir_outputs_mc.ok"),
        log:
            "{folder}/log/publishdir_outputs/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        script:
            "../scripts/utils/publishdir.py"
