# Pipeline input rules

rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        ancient(config["output_location"] + "config/config_df.tsv"),
        bam = config["input_bam_location"]
    output:
        config["output_location"] + "config/exclude_file"
    params:
        chroms = config["chromosomes"]
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/generate_exclude_file.py"


