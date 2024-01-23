# if config["multistep_normalisation"] == False or config["ashleys_pipeline"] == False:


rule mergeBams:
    input:
        check=remove_unselected_fct,
        bam=selected_input_bam,
        bai=selected_input_bai,
        labels="{folder}/{sample}/cell_selection/labels.tsv",
    output:
        temp("{folder}/{sample}/merged_bam/merged.raw.bam"),
    log:
        "{folder}/log/mergeBams/{sample}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time=60,
    threads: 10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools merge -@ {threads} {output} {input.bam} 2>&1 > {log}"


rule mergeSortBams:
    input:
        "{folder}/{sample}/merged_bam/merged.raw.bam",
    output:
        temp("{folder}/{sample}/merged_bam/merged.bam"),
    log:
        "{folder}/log/mergeBams/{sample}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time=60,
    threads: 10
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"


rule index_merged_bam:
    input:
        "{folder}/{sample}/merged_bam/merged.bam",
    output:
        temp("{folder}/{sample}/merged_bam/merged.bam.bai"),
    log:
        "{folder}/log/merged_bam/{sample}/merged.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "samtools index {input} > {log} 2>&1"


rule regenotype_SNVs:
    input:
        bam="{folder}/{sample}/merged_bam/merged.bam",
        bai="{folder}/{sample}/merged_bam/merged.bam.bai",
        sites=config["references_data"][config["reference"]]["snv_sites_to_genotype"],
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(
            fasta=config["references_data"][config["reference"]]["reference_fasta"]
        ),
    output:
        vcf="{folder}/{sample}/snv_genotyping/{chrom,chr[0-9A-Z]+}.vcf",
    log:
        "{folder}/log/snv_genotyping/{sample}/{chrom}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
        time=600,
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
        bam="{folder}/{sample}/merged_bam/merged.bam",
        bai="{folder}/{sample}/merged_bam/merged.bam.bai",
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(
            fasta=config["references_data"][config["reference"]]["reference_fasta"]
        ),
        ploidy="{folder}/{sample}/ploidy/ploidy_bcftools.txt",
    output:
        vcf="{folder}/{sample}/snv_calls/{chrom,chr[0-9A-Z]+}.vcf",
    log:
        "{folder}/log/snv_calls/{sample}/{chrom,chr[0-9A-Z]+}.vcf",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time=600,
    shell:
        """
        bcftools mpileup -r {wildcards.chrom} -f {input.fasta} {input.bam} \
        | bcftools call -mv --ploidy-file {input.ploidy} | bcftools view --genotype het --types snps > {output} 2> {log}
        """


checkpoint check_SNVs_nb:
    input:
        lambda wc: expand(
            "{folder}/{sample}/snv_calls/{chrom}.vcf",
            chrom=config["chromosomes"],
            sample=wc.sample,
            folder=wc.folder,
        ),
    output:
        summary_snp_nb="{folder}/{sample}/snv_calls/check_SNVs_nb.txt",
    log:
        "{folder}/log/check_SNVs_nb/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        # header
        echo "chrom\tSNP_nb" >> {output.summary_snp_nb}
        # for each chromosome, check the number of SNVs
        for chrom_input in {input};
        do
            chrom=$(basename $chrom_input .vcf);
            nb=$(grep -v "#" $chrom_input | wc -l);
            echo "$chrom\t$nb" >> {output.summary_snp_nb};
        done
        """
