
rule prepare_haplotag_input_bam:
    """
    Materialise the BAM consumed by haplotag_bams with an @RG SM tag matching the
    sample. whatshap matches VCF sample name <-> BAM SM tag; pooled libraries carry
    the pool name (e.g. embl_9i) instead of the donor (e.g. HG01412).
      SM correct                     -> hardlink (no disk cost)
      SM wrong + fix_sm_tag: True    -> samtools reheader rewrites SM to {sample}
      SM wrong + fix_sm_tag: False   -> fail loudly, naming both values
    """
    input:
        bam=lambda wc: expand(
            "{input_folder}/{{sample}}/selected/{{cell}}.sort.mdup.bam",
            input_folder=config["data_location"],
        ),
        bai=lambda wc: expand(
            "{input_folder}/{{sample}}/selected/{{cell}}.sort.mdup.bam.bai",
            input_folder=config["data_location"],
        ),
    output:
        bam=temp("{folder}/{sample}/haplotag/input/{cell}.bam"),
        bai=temp("{folder}/{sample}/haplotag/input/{cell}.bam.bai"),
    log:
        "{folder}/log/prepare_haplotag_input_bam/{sample}/{cell}.log",
    container:
        None
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    envmodules:
        "SAMtools/1.21-GCC-13.3.0",
    threads: 1
    resources:
        mem_mb=2000,
    params:
        fix=config.get("fix_sm_tag", False),
    shell:
        r"""
        exec > {log} 2>&1
        set -o pipefail
        sm=$(samtools view -H {input.bam} | grep -m1 '^@RG' | tr '\t' '\n' | sed -n 's/^SM://p' | head -1)
        echo "cell={wildcards.cell} sample={wildcards.sample} SM=$sm fix={params.fix}"
        if [ "$sm" = "{wildcards.sample}" ]; then
            if cp -l {input.bam} {output.bam} 2>/dev/null && cp -l {input.bai} {output.bai} 2>/dev/null; then
                echo "SM matches sample -> hardlink"
            else
                echo "SM matches sample -> cross-device, symlink"
                rm -f {output.bam} {output.bai}
                ln -s "$(readlink -f {input.bam})" {output.bam}
                ln -s "$(readlink -f {input.bai})" {output.bai}
            fi
        elif [ "{params.fix}" = "True" ]; then
            echo "SM mismatch -> reheader SM:$sm to SM:{wildcards.sample}"
            samtools reheader -c "sed 's/\tSM:[^\t]*/\tSM:{wildcards.sample}/'" {input.bam} > {output.bam}
            samtools index {output.bam}
        else
            echo "ERROR: BAM carries SM:$sm but sample is {wildcards.sample}."
            echo "       whatshap would fail with 'No common samples between VCF and BAM'."
            echo "       Set fix_sm_tag: True (or --config fix_sm_tag=True) to rewrite the tag."
            exit 1
        fi

        # A hardlink shares its inode with the input and a symlink resolves to it, so the
        # output inherits the input .bam's mtime - which is older than the .bai written a
        # moment later, and snakemake fails any output older than an input. Stamping every
        # path in ONE touch call gives them all the same mtime, so none is ever older.
        touch {output.bam} {output.bai}
        """


rule haplotag_bams:
    input:
        vcf="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz",
        tbi="{folder}/{sample}/strandphaser/phased-snvs/{sample}.vcf.gz.tbi",
        bam="{folder}/{sample}/haplotag/input/{cell}.bam",
        bai="{folder}/{sample}/haplotag/input/{cell}.bam.bai",
        # check="{folder}/{sample}/config/remove_unselected_bam.ok",
        # check=remove_unselected_fct,
        # bam=selected_input_bam,
        # bam=lambda wc: expand(
        #     "{input_folder}/{{sample}}/bam/{{cell}}.sort.mdup.bam",
        #     input_folder=config["data_location"],
        # ),
        fasta=get_reference_fasta(),
        fasta_index=f"{get_reference_fasta()}.fai",
    output:
        "{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
    log:
        "{folder}/log/haplotag_bams/{sample}/{cell}.log",
    group:
        "haplotagging_per_cell"
    params:
        ref=config["reference"],
    resources:
        mem_mb=get_mem_mb_haplotag_group,
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "whatshap haplotag --skip-missing-contigs -o {output} -r {input.fasta} {input.vcf} {input.bam} > {log} 2>{log}  "


rule create_haplotag_segment_bed:
    localrule: True
    group:
        "text_processing_operations"
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
        mem_mb=get_mem_mb_create_haplotag_table,
    script:
        "../scripts/haplotagging_scripts/haplotagTable.snakemake.R"


rule merge_haplotag_tables:
    localrule: True
    input:
        # tsvs=lambda wc: [
        #     "{}/{}/haplotag/table/by-cell/{}.tsv".format(
        #         config["data_location"], wc.sample, cell
        #     )
        #     for cell in bam_per_sample[wc.sample]
        # ],
        tsvs=aggregate_cells_haplotag_tables,
    output:
        tsv="{folder}/{sample}/haplotag/table/haplotag_counts_merged.tsv",
    log:
        "{folder}/log/haplotag/table/{sample}/haplotag_counts_merged.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        # Convert Snakemake's space-separated string to bash array
        TSVS=({input.tsvs})

        if [ ${{#TSVS[@]}} -eq 0 ]; then
            echo "WARNING: No cells passed QC filter - creating empty haplotag table" >&2
            # Create header-only output for downstream compatibility
            echo -e "cell\\tchrom\\tstart\\tend\\tcrick.count\\twatson.count\\tcrick.H1\\tcrick.H2\\twatson.H1\\twatson.H2" > {output.tsv}
        else
            (head -n1 ${{TSVS[0]}} && tail -q -n +2 {input.tsvs}) > {output.tsv}
        fi
        """
