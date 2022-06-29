# Pipeline input rules

rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        # ancient(config["output_location"] + "config/config_df.tsv"),
        ancient("config/samples.tsv"),
        bam = config["input_bam_location"]
    output:
        "{output}/config_output/exclude_file"
    log:
        "{output}/log/config_output/exclude_file.log"
    params:
        chroms = config["chromosomes"]
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/generate_exclude_file.py"


