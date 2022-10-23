
rule haplotag_bams:
    input:
        vcf="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz",
        tbi="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
        bam=lambda wc: expand(
            "{input_folder}/{{sample}}/bam/{{cell}}.sort.mdup.bam",
            input_folder=config["data_location"],
        ),
        fasta=config["references_data"][config["reference"]]["reference_fasta"],
        fasta_index="{fasta}.fai".format(
            fasta=config["references_data"][config["reference"]]["reference_fasta"]
        ),
    output:
        "{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
    log:
        "{folder}/log/haplotag_bams/{sample}/{cell}.log",
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
        segments="{folder}/{sample}/segmentation/Selection_jointseg.txt",
    output:
        bed="{folder}/{sample}/haplotag/bed/{sample}.bed",
    log:
        "{folder}/log/haplotag/bed/{sample}.log",
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
        bam="{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
        bai="{folder}/{sample}/haplotag/bam/{cell}.bam.htg.bai",
        bed="{folder}/{sample}/haplotag/bed/{sample}.bed",
        paired_end="{folder}/{sample}/config/single_paired_end_detection.txt",
    output:
        tsv="{folder}/{sample}/haplotag/table/by-cell/{cell}.tsv",
    log:
        "{folder}/log/create_haplotag_table/{sample}.{cell}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/haplotagging_scripts/haplotagTable.snakemake.R"


rule merge_haplotag_tables:
    input:
        tsvs=lambda wc: [
            "{}/{}/haplotag/table/by-cell/{}.tsv".format(
                config["data_location"], wc.sample, cell
            )
            for cell in bam_per_sample[wc.sample]
        ],
    output:
        tsv="{folder}/{sample}/haplotag/table/haplotag_counts_merged.tsv",
    log:
        "{folder}/log/haplotag/table/{sample}/haplotag_counts_merged.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        "(head -n1 {input.tsvs[0]} && tail -q -n +2 {input.tsvs}) > {output.tsv}"
