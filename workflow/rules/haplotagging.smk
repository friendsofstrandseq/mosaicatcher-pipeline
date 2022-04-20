import pandas as pd
config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# print(config_df)
bam_per_sample = config_df.loc[config_df["all/selected"] == "selected"].groupby("Sample")["File"].apply(list).to_dict()


################################################################################
# Haplotagging                                                                 #
################################################################################



rule haplotag_bams:
    input:
        vcf = config["output_location"] + "strandphaser/phased-snvs/{sample}.vcf.gz",
        tbi = config["output_location"] + "strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
        bam = config["input_bam_location"] + "{sample}/selected/{cell}.bam",
        bai = config["input_bam_location"] + "{sample}/selected/{cell}.bam.bai"
    output:
        config["output_location"] + "haplotag/bam/{sample}/{cell}.bam",
    log:
        config["output_location"] + "log/haplotag_bams/{sample}/{cell}.log"
    params:
        ref = config["reference"]
    conda: 
        "../envs/mc_bioinfo_tools.yaml"
    shell:
    # DOCME : --skip-missing-contigs option to remove unused chroms
        # "whatshap haplotag -o {output} -r {params.ref} {input.vcf} {input.bam} > {log} 2>{log}  " 
        "whatshap haplotag --skip-missing-contigs -o {output} -r {params.ref} {input.vcf} {input.bam} > {log} 2>{log}  " 

rule create_haplotag_segment_bed:
    input:
        segments = config["output_location"] + "segmentation/{sample}/Selection_jointseg.txt",
    output:
        bed = config["output_location"] + "haplotag/bed/{sample}.bed",
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v s={config[window]} -f scripts/haplotagging_scripts/create_haplotag_segment_bed.awk {input.segments} > {output.bed}
        """

rule create_haplotag_table:
    input:
        bam = config["output_location"] + "haplotag/bam/{sample}/{cell}.bam",
        bai = config["output_location"] + "haplotag/bam/{sample}/{cell}.bam.bai",
        bed = config["output_location"] + "haplotag/bed/{sample}.bed"
    output:
        tsv = config["output_location"] + "haplotag/table/{sample}/by-cell/{cell}.tsv"
    log:
        config["output_location"] + "log/create_haplotag_table/{sample}.{cell}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/haplotagging_scripts/haplotagTable.snakemake.R"

rule merge_haplotag_tables:
    input:
        tsvs = lambda wc: [config["output_location"] + "haplotag/table/{}/by-cell/{}.tsv".format(wc.sample, cell) for cell in bam_per_sample[wc.sample]],
    output:
        tsv = config["output_location"] + "haplotag/table/{sample}/haplotag_counts_merged.tsv"
    shell:
        "(head -n1 {input.tsvs[0]} && tail -q -n +2 {input.tsvs}) > {output.tsv}"


