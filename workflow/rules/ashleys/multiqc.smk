
rule ashleys_fastqc:
    input:
        "{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
    output:
        html=report(
            "{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.html",
            category="FastQC",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Cell": "{cell}", "Pair": "{pair}"},
        ),
        zip="{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.zip",
    log:
        "{folder}/log/fastqc/{sample}/{cell}_{pair}.log",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "FastQC/0.11.9-Java-11",
    shell:
        """
        tempdir=$(mktemp -d)
        fastqc --quiet -t {threads} --outdir "$tempdir" {input} > {log} 2>&1
        mv "$tempdir/{wildcards.cell}.{wildcards.pair}_fastqc.html" {output.html}
        mv "$tempdir/{wildcards.cell}.{wildcards.pair}_fastqc.zip" {output.zip}
        rm -rf "$tempdir"
        """


rule ashleys_fastqc_aggregate:
    localrule: True
    input:
        lambda wc: expand(
            "{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.html",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
            pair=[1, 2],
        ),
    output:
        touch("{folder}/{sample}/multiqc/fastqc/config/fastqc_output_touch.ok"),


rule ashleys_samtools_idxstats:
    group: "qc_statistics_operations"
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_idxstats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_idxstats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
    shell:
        "samtools idxstats {input} > {output}"


rule ashleys_samtools_idxstats_aggr:
    localrule: True
    input:
        lambda wc: expand(
            "{folder}/{sample}/multiqc/samtools_idxstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok"
        ),


rule ashleys_samtools_flagstats:
    group: "qc_statistics_operations"
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_flagstats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_flagstats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
    shell:
        "samtools flagstats {input} > {output}"


rule ashleys_samtools_flagstats_aggr:
    localrule: True
    input:
        lambda wc: expand(
            "{folder}/{sample}/multiqc/samtools_flagstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok"
        ),


rule ashleys_samtools_stats:
    group: "qc_statistics_operations"
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_stats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_stats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
    shell:
        "samtools stats {input} > {output}"


rule ashleys_samtools_stats_aggr:
    localrule: True
    input:
        lambda wc: expand(
            "{folder}/{sample}/multiqc/samtools_stats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_stats/config/samtools_stats_aggr_touch.ok"
        ),


rule ashleys_multiqc:
    input:
        fastqc="{folder}/{sample}/multiqc/fastqc/config/fastqc_output_touch.ok",
        samtools_idxstats="{folder}/{sample}/multiqc/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok",
        samtools_stats="{folder}/{sample}/multiqc/samtools_stats/config/samtools_stats_aggr_touch.ok",
        samtools_flagstats="{folder}/{sample}/multiqc/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok",
    output:
        report="{folder}/{sample}/multiqc/multiqc_report/multiqc_report.html",
        outdir=report(
            directory("{folder}/{sample}/multiqc/multiqc_report"),
            htmlindex="multiqc_report.html",
            category="MultiQC",
            subcategory="{sample}",
        ),
    log:
        "{folder}/{sample}/log/multiqc/{sample}.log",
    resources:
        mem_mb=get_mem_mb_heavy,
    params:
        multiqc_input=lambda wc, input: "{abs_path}".format(
            abs_path=config["abs_path"]
        ).join(input.fastqc.split("/")[:-3]),
    conda:
        "../../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "MultiQC/1.30-foss-2024a",
    shell:
        "multiqc {params.multiqc_input} --outdir {output.outdir}"
