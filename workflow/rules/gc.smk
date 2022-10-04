if config["GC_analysis"] is True:

    rule alfred:
        input:
            merged_bam="{output_folder}/merged_bam/{sample}/merged.bam",
            merged_bam_bai="{output_folder}/merged_bam/{sample}/merged.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{output_folder}/alfred/{sample}.json.gz",
            alfred_tsv="{output_folder}/alfred/{sample}.tsv.gz",
        log:
            "{output_folder}/log/alfred/{sample}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.merged_bam}
            """

    rule alfred_table:
        input:
            "{output_folder}/alfred/{sample}.tsv.gz",
        output:
            "{output_folder}/alfred/{sample}.table",
        log:
            "{output_folder}/log/alfred_table/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            zcat {input} | grep "^GC" > {output}
            """

    rule VST_correction:
        input:
            counts="{output_folder}/counts/{sample}/{sample}.txt.filter.gz",
        output:
            counts_vst="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.txt.gz",
        log:
            "{output_folder}/log/VST_correction/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/variance_stabilizing_transformation.R"

    rule GC_correction:
        input:
            counts_vst="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.txt.gz",
        output:
            counts_vst_gc="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.GC.txt.gz",
        log:
            "{output_folder}/log/GC_correction/{sample}.log",
        params:
            gc_matrix="workflow/data/GC/GC_matrix_200000.txt",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/GC_correction.R"

    rule counts_scaling:
        input:
            counts_vst_gc="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.GC.txt.gz",
        output:
            counts_vst_gc_scaled="{output_folder}/counts/GC_correction/{sample}/{sample}.VST.GC.scaled.txt.gz",
        log:
            "{output_folder}/log/counts_scaling/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/counts_scaling.R"
