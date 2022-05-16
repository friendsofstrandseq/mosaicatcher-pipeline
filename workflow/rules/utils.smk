
################################################################################
# UTILS                                                                        #
################################################################################



# rule sort_bam:
#     input:
#         config["output_location"] + "snv_calls/{sample}/merged.unsorted.bam"
#     output:
#         config["output_location"] + "snv_calls/{sample}/merged.bam"
#     shell:
#         config["samtools"] + " sort {input} -o {output}"


rule index_bam:
    """
    rule fct:
    input:
    output:
    """
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    # log:
    #     "{file}.bam.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        # config["samtools"] + " index {input} 2> {log}"
        "samtools" + " index {input}"


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
    # log:
    #     "log/compress_vcf/{file}.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        # "(cat {input.vcf} | bgzip > {output.vcf}) > {log} 2>&1"
        "bgzip {input.vcf}"

rule index_vcf:
    """
    rule fct:
    input:
    output:
    """
    input:
        vcf="{file}.vcf.gz",
    output:
        tbi="{file}.vcf.gz.tbi",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "tabix -p vcf {input.vcf}"

