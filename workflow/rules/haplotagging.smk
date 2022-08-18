# from workflow.scripts.utils.utils import get_mem_mb

# import pandas as pd
# config_df = pd.read_csv("config/config_df.tsv", sep="\t")
# bam_per_sample = config_df.loc[config_df["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()
# bam_per_sample = df_config_files.loc[df_config_files["Selected"] == True].groupby("Sample")["File"].apply(list).to_dict()

################################################################################
# Haplotagging                                                                 #
################################################################################

# DOCME : --skip-missing-contigs option to remove unused chroms
# "whatshap haplotag -o {output_folder} -r {params.ref} {input.vcf} {input.bam} > {log} 2>{log}  "
# bai = config["input_bam_location"] + "{sample}/selected/{cell}.bam.bai"


rule haplotag_bams:
    input:
        vcf="{output_folder}/strandphaser/phased-snvs/{sample}.vcf.gz",
        tbi="{output_folder}/strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
        bam=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam",
            zip,
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)]
            if wc.sample in bam_per_sample_local
            else "FOOBAR",
        ),
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(fasta=config["references_data"][config["reference"]]["reference_fasta"]),
    output:
        "{output_folder}/haplotag/bam/{sample}/{cell}.bam",
    log:
        "{output_folder}/log/haplotag_bams/{sample}/{cell}.log",
    params:
        ref=config["reference"],
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "whatshap haplotag --skip-missing-contigs -o {output} -r {input.fasta} {input.vcf} {input.bam} > {log} 2>{log}  "


rule create_haplotag_segment_bed:
    input:
        segments="{output_folder}/segmentation/{sample}/Selection_jointseg.txt",
    output:
        bed="{output_folder}/haplotag/bed/{sample}.bed",
    log:
        "{output_folder}/log/haplotag/bed/{sample}.log",
    params:
        window=config["window"],
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
        awk -v s={params.window} -f workflow/scripts/haplotagging_scripts/create_haplotag_segment_bed.awk {input.segments} > {output.bed}
        """


rule create_haplotag_table:
    input:
        bam="{output_folder}/haplotag/bam/{sample}/{cell}.bam",
        bai="{output_folder}/haplotag/bam/{sample}/{cell}.bam.bai",
        bed="{output_folder}/haplotag/bed/{sample}.bed",
        paired_end="{output_folder}/config/{sample}/single_paired_end_detection.txt",
    output:
        tsv="{output_folder}/haplotag/table/{sample}/by-cell/{cell}.tsv",
    log:
        "{output_folder}/log/create_haplotag_table/{sample}.{cell}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/haplotagging_scripts/haplotagTable.snakemake.R"


rule merge_haplotag_tables:
    input:
        tsvs=lambda wc: [
            "{}/haplotag/table/{}/by-cell/{}.tsv".format(
                config["output_location"], wc.sample, cell
            )
            for cell in bam_per_sample[wc.sample]
        ],
    output:
        tsv="{output_folder}/haplotag/table/{sample}/haplotag_counts_merged.tsv",
    log:
        "{output_folder}/log/haplotag/table/{sample}/haplotag_counts_merged.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "(head -n1 {input.tsvs[0]} && tail -q -n +2 {input.tsvs}) > {output.tsv}"
