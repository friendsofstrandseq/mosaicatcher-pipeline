## Rules to perform GC analysis & correction on Strand-Seq libraries
## ---------------------------------------------------------------
## fastqc: QC analysis of FASTQ files
## bwa_index: index fasta file using the BW transformation to perform next reads mapping
## bwa_strandseq_to_reference_alignment: mapping against FASTA reference (based on reference selected (hg19/hg38/T2T))
## ashleys_samtools_sort_bam: sorting bam files
## ashleys_mark_duplicates: mark duplicates in bam files
## *samtools index: index bam files (only if not loaded as a module into mosaicatcher_pipeline)
## generate_features/predict: features creation & prediction using ashleys-qc ML method to detect high/low quality libraries
## notebook_hand_selection: fire a jupyter notebook that allow hand selection of low quality cells based on QC plots


if config["genecore"] is True and config["genecore_date_folder"]:
    if config["mosaicatcher_pipeline"] is False:

        localrules:
            ashleys_genecore_symlink,
            ashleys_symlink_bam_ashleys,

    rule ashleys_genecore_symlink:
        input:
            lambda wc: df_config_files.loc[
                (df_config_files["Sample"] == wc.sample)
                & (df_config_files["File"] == "{}.{}".format(wc.cell, wc.pair))
            ]["Genecore_path"]
            .unique()
            .tolist(),
        output:
            "{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
        log:
            "{folder}/log/genecore_symlink/{sample}/{cell}_{pair}.log",
        shell:
            "ln -f -s {input} {output}"

    ruleorder: ashleys_genecore_symlink > ashleys_bwa_strandseq_to_reference_alignment


localrules:
    ashleys_symlink_bam_ashleys,


rule ashleys_bwa_index:
    input:
        ancient(config["references_data"][config["reference"]]["reference_fasta"]),
    output:
        idx=multiext(
            config["references_data"][config["reference"]]["reference_fasta"],
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{}.log".format(
            config["references_data"][config["reference"]]["reference_fasta"]
        ),
    params:
        algorithm="bwtsw",
    threads: 16
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/bwa/index"


ruleorder: ashleys_bwa_strandseq_to_reference_alignment > ashleys_samtools_sort_bam > ashleys_mark_duplicates


if config["paired_end"] is True:

    rule ashleys_bwa_strandseq_to_reference_alignment:
        input:
            mate1="{folder}/{sample}/fastq/{cell}.1.fastq.gz",
            mate2="{folder}/{sample}/fastq/{cell}.2.fastq.gz",
            ref="{ref}".format(
                ref=config["references_data"][config["reference"]]["reference_fasta"]
            ),
            ref_index=multiext(
                config["references_data"][config["reference"]]["reference_fasta"],
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            ),
        output:
            bam=temp("{folder}/{sample}/bam/{cell}.bam.raw"),
        log:
            bwa="{folder}/{sample}/log/{cell}.bwa.log",
            samtools="{folder}/{sample}/log/{cell}.samtools.log",
        threads: 6
        params:
            idx_prefix=lambda wildcards, input: input.ref_index[0].rsplit(".", 1)[0],
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        conda:
            "../../envs/mc_base.yaml"
        shell:
            "bwa mem -t {threads}"
            ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
            " -v 2 {input.ref} {input.mate1} {input.mate2} 2> {log.bwa} | "
            " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"

else:

    rule ashleys_bwa_strandseq_to_reference_alignment:
        input:
            mate1="{folder}/{sample}/fastq/{cell}.1.fastq.gz",
            ref="{ref}".format(
                ref=config["references_data"][config["reference"]]["reference_fasta"]
            ),
            ref_index=multiext(
                config["references_data"][config["reference"]]["reference_fasta"],
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            ),
        output:
            bam=temp("{folder}/{sample}/bam/{cell}.bam.raw"),
        log:
            bwa="{folder}/{sample}/log/{cell}.bwa.log",
            samtools="{folder}/{sample}/log/{cell}.samtools.log",
        threads: 6
        params:
            idx_prefix=lambda wildcards, input: input.ref_index[0].rsplit(".", 1)[0],
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        conda:
            "../../envs/mc_base.yaml"
        shell:
            "bwa mem -t {threads}"
            ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
            " -v 2 {input.ref} {input.mate1} 2> {log.bwa} | "
            " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"


rule ashleys_samtools_sort_bam:
    input:
        "{folder}/{sample}/bam/{cell}.bam.raw",
    output:
        temp("{folder}/{sample}/bam/{cell}.bam.sort"),
    log:
        "{folder}/{sample}/log/samtools_sort/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    conda:
        "../../envs/mc_base.yaml"
    shell:
        "samtools sort -O BAM -o {output} {input} 2>&1 > {log}"


rule ashleys_mark_duplicates:
    input:
        bam="{folder}/{sample}/bam/{cell}.bam.sort",
    output:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    log:
        "{folder}/{sample}/log/markdup/{cell}.log",
    conda:
        "../../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="01:00:00",
    shell:
        "sambamba markdup {input.bam} {output} 2>&1 > {log}"


rule ashleys_symlink_bam_ashleys:
    input:
        bam="{folder}/{sample}/bam/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    output:
        bam="{folder}/{sample}/bam_ashleys/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam_ashleys/{cell}.sort.mdup.bam.bai",
    log:
        "{folder}/log/symlink_selected_bam/{sample}/{cell}.log",
    conda:
        "../../envs/mc_base.yaml"
    script:
        "../../scripts/ashleys/utils/symlink_selected_bam.py"
