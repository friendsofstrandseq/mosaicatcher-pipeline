

rule generate_exclude_file_for_mosaic_count:
    """
    rule fct: 
    input:
    output:
    """
    input:
        # ancient("config/samples.tsv"),
        bam=lambda wc: expand(
            "{input_folder}/{sample}/bam/{cell}.sort.mdup.bam",
            input_folder=config["data_location"],
            sample=wc.sample,
            cell=bam_per_sample_local[str(wc.sample)],
        ),
    output:
        "{folder}/{sample}/config/chroms_to_exclude.txt",
    log:
        "{folder}/log/config_output/{sample}/exclude_file.log",
    params:
        chroms=config["chromosomes"],
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/generate_exclude_file.py"
