# rule prepare_breakpointr:
#     input:
#         seg_initial_str_state="{folder}/{sample}/segmentation/Selection_initial_strand_state",
#         single_paired_end_detect="{folder}/{sample}/config/single_paired_end_detection.txt",
#     output:
#         # "{folder}/{sample}/strandphaser/StrandPhaseR.{chrom}.config",
#     log:
#         # "{folder}/log/strandphaser/{sample}/StrandPhaseR.{chrom}.log",
#     conda:
#         "../envs/mc_base.yaml"
#     script:
#         "../scripts/breakpointr_scripts/prepare_breakpointr.py"


rule run_breakpointr:
    input:
        bsgenome="workflow/data/ref_genomes/config/BSgenome_{}.ok".format(
            config["reference"]
        ),
        bams=selected_input_bam,
    output:
        output_folder=directory("{folder}/{sample}/breakpointR"),
        output_config="{folder}/{sample}/breakpointR/breakpointR.config",
    log:
        "{folder}/log/run_breakpointr/{sample}.log",
    conda:
        "../envs/breakpointr.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time=600,
    params:
        input_bam=lambda wc: "{}/{}/selected".format(config["data_location"], wc.sample),
    shell:
        """
        Rscript workflow/scripts/breakpointr_scripts/breakpointR_pipeline.R \
                {params.input_bam} \
                {output.output_folder}
        """
