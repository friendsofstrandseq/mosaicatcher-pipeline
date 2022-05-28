# Pipeline input rules

rule check_bam_input:
    """
    rule fct:
    input:
    output:
    """
    input:
        config["input_bam_location"]
    output:
        config["output_location"] + "config/config_df.tsv"
    conda:
        "../envs/mc_base.yaml"
    params:
        check_sm_tag = config["check_sm_tag"]
    script:
        "../scripts/utils/handle_input.py"