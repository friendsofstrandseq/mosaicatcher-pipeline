import pandas as pd
config_df = pd.read_csv(config["output_location"] + "config/config_df.tsv", sep="\t")
allbams_per_sample = df_config_files.loc[df_config_files["all/selected"] == "all"].groupby("Sample")["File"].apply(list).to_dict()


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
    threads:
        4
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        # FIXME : Samtools 1.10 from Conda env not working ; 1.9 from Seneca working > change it into conda env yml file
        # config["samtools"] + " merge -@ {threads} {output} {input} 2>&1 > {log}"
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
        # bcftools = config["bcftools"]
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
    # CHECKME : Samtools / BCFtools / freebayes path definition through conda env
    # CHECKME : interest of using -r parameters for freebayes => split by chroms
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