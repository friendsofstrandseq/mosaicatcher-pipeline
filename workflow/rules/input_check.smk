# Pipeline input rules



rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        # ancient("config/samples.tsv"),
        bam=lambda wc: expand(
            "{input_folder}/{sample}/all/{cell}.sort.mdup.bam",
            input_folder=config["input_bam_location"],
            sample=samples,
            cell=bam_per_sample_local[str(wc.sample)]
            if wc.sample in bam_per_sample_local
            else "FOOBAR",
        ),
    output:
        "{output_folder}/config_output/{sample}/exclude_file",
    log:
        "{output_folder}/log/config_output/{sample}/exclude_file.log",
    params:
        chroms=config["chromosomes"],
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/generate_exclude_file.py"
