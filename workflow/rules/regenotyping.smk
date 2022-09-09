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
        lambda wc: expand(
            "{input_folder}/{sample}/all/{bam}.sort.mdup.bam",
            input_folder=config["input_bam_location"],
            sample=wc.sample,
            bam=allbams_per_sample[wc.sample],
        ),
    output:
        "{output_folder}/merged_bam/{sample}/merged.raw.bam",
    log:
        "{output_folder}/log/mergeBams/{sample}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="01:00:00",
    threads: 10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"


rule mergeSortBams:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{output_folder}/merged_bam/{sample}/merged.raw.bam",
    output:
        "{output_folder}/merged_bam/{sample}/merged.bam",
    log:
        "{output_folder}/log/mergeBams/{sample}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="01:00:00",
    threads: 10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"


rule regenotype_SNVs:
    """
    rule fct:
    input:
    output:
    """
    input:
        bam="{output_folder}/merged_bam/{sample}/merged.bam",
        bai="{output_folder}/merged_bam/{sample}/merged.bam.bai",
        # sites=config["snv_sites_to_genotype"],
        sites=config["references_data"][config["reference"]]["snv_sites_to_genotype"],
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(
            fasta=config["references_data"][config["reference"]]["reference_fasta"]
        ),
    output:
        vcf="{output_folder}/snv_genotyping/{sample}/{chrom,chr[0-9A-Z]+}.vcf",
    log:
        "{output_folder}/log/snv_genotyping/{sample}/{chrom}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        """
        (freebayes \
            -f {input.fasta} \
            -r {wildcards.chrom} \
            -@ {input.sites} \
            --only-use-input-alleles {input.bam} \
            --genotype-qualities \
        | bcftools view \
            --exclude-uncalled \
            --types snps \
            --genotype het \
            --include "QUAL>=10" \
        > {output.vcf}) 2> {log}
        """


rule call_SNVs_bcftools_chrom:
    input:
        bam="{output_folder}/merged_bam/{sample}/merged.bam",
        bai="{output_folder}/merged_bam/{sample}/merged.bam.bai",
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(
            fasta=config["references_data"][config["reference"]]["reference_fasta"]
        ),
        ploidy="{output_folder}/ploidy/{sample}/ploidy_bcftools.txt",
    output:
        vcf="{output_folder}/snv_calls/{sample}/{chrom,chr[0-9A-Z]+}.vcf",
    log:
        "{output_folder}/log/snv_calls/{sample}/{chrom,chr[0-9A-Z]+}.vcf",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    shell:
        """
        bcftools mpileup -r {wildcards.chrom} -f {input.fasta} {input.bam} \
        | bcftools call -mv --ploidy-file {input.ploidy} | bcftools view --genotype het --types snps > {output} 2> {log}
        """
