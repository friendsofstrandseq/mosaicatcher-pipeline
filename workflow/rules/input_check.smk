# Pipeline input rules


rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        config["output_location"] + "config/config_df.tsv",
        bam = config["input_bam_location"]
    output:
        config["output_location"] + "config/exclude_file"
    params:
        chroms = config["chromosomes"]
    conda:
        "../envs/mc_base.yaml"
    # shell:
    #     "python scripts/utils/generate_exclude_file.py {input} {output} {params.chroms}"
    script:
        "../scripts/utils/generate_exclude_file.py"


