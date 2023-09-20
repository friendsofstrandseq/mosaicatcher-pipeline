
rule fastqc_debug:
    input:
        fastq="{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
        fastqc_check="{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.html",
        labels_strandscape="{folder}/{sample}/cell_selection/labels_strandscape.tsv",
    output:
        html=report(
            "{folder}/{sample}/debug/mosaicatcher_fastqc/{cell}_{pair}_fastqc.html",
            category="FastQC",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Cell": "{cell}", "Pair": "{pair}"},
        ),
        zip="{folder}/{sample}/debug/mosaicatcher_fastqc/{cell}_{pair}_fastqc.zip",
    log:
        "{folder}/log/fastqc_debug/{sample}/{cell}_{pair}.log",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/dev/debug.yaml"
    params:
        outdir = lambda wc, output: "/".join(output.zip.split("/")[:-1])
    shell:
        "fastqc --outdir {params.outdir} --quiet {input.fastq}"