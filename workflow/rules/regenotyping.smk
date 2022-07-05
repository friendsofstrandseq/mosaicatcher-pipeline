from workflow.scripts.utils.utils import get_mem_mb 

import pandas as pd
config_df = pd.read_csv(config["output_location"] + "config/config_df.tsv", sep="\t")
allbams_per_sample = config_df.groupby("Sample")["File"].apply(list).to_dict()


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
        lambda wc: expand(config["input_bam_location"] + wc.sample + "/all/{bam}.bam", bam = allbams_per_sample[wc.sample]) if wc.sample in allbams_per_sample else "FOOBAR",
    output:
        config["output_location"] + "merged_bam/{sample}/merged.bam"
    log:
        config["output_location"] + "log/mergeBams/{sample}.log"
    resources:
        mem_mb = get_mem_mb,
        time = "01:00:00",
    threads:
        10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools" + " merge -@ {threads} {output} {input} 2>&1 > {log}"

rule regenotype_SNVs:
    """
    rule fct:
    input:
    output:
    """
    input:
        bam   = config["output_location"] + "merged_bam/{sample}/merged.bam",
        bai   = config["output_location"] + "merged_bam/{sample}/merged.bam.bai",
        sites = config["snv_sites_to_genotype"],
    output:
        vcf = config["output_location"] + "snv_genotyping/{sample}/{chrom,chr[0-9A-Z]+}.vcf"
    log:
        config["output_location"] + "log/snv_genotyping/{sample}/{chrom}.log"
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
            --types snps \
            --genotype het \
            --include "QUAL>=10" \
        > {output.vcf}) 2> {log}
        """

