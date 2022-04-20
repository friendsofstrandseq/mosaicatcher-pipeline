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
        "config/config_df.tsv"
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/handle_input.py"