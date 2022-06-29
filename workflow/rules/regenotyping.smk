
################################################################################
# REGENOTYPE SNV                                                               #
################################################################################

rule mergeBams:
    """
    rule fct:
    input:
    output:
    """
    input:
        lambda wc: expand("{input_folder}/{sample}/all/{bam}.bam", input_folder=config['input_bam_location'], sample=samples, bam = allbams_per_sample[wc.sample]) if wc.sample in allbams_per_sample else "FOOBAR",
    output:
        "{output}/merged_bam/{sample}/merged.bam"
    log:
        "{output}/log/mergeBams/{sample}.log"
    resources:
        mem_mb = get_mem_mb,
        time = "01:00:00",
    threads:
        10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"

rule regenotype_SNVs:
    """
    rule fct:
    input:
    output:
    """
    input:
        bam   = "{output}/merged_bam/{sample}/merged.bam",
        bai   = "{output}/merged_bam/{sample}/merged.bam.bai",
        sites = config["snv_sites_to_genotype"],
    output:
        vcf = "{output}/snv_genotyping/{sample}/{chrom,chr[0-9A-Z]+}.vcf"
    log:
        "{output}/log/snv_genotyping/{sample}/{chrom}.log"
    params:
        fa = config["reference"],
    resources:
        mem_mb = get_mem_mb,
        time = "10:00:00",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        """
        (freebayes \
            -f {params.fa} \
            -r {wildcards.chrom} \
            -@ {input.sites} \
            --only-use-input-alleles {input.bam} \
            --genotype-qualities \
        | bcftools view \
            --exclude-uncalled \
            --genotype het \
            --types snps \
            --include "QUAL>=10" \
        > {output.vcf}) 2> {log}
        """

