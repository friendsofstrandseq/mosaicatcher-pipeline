
rule install_rlib_strandphaser:
    output:
        check=touch(
            expand(
                "{}/strandphaser/R_setup/strandphaser_version-{}.ok".format(
                    config["output_location"], config["git_commit_strandphaser"]
                )
            )
        ),
    log:
        expand(
            "{}/log/strandphaser/R_setup/strandphaser_version-{}.ok".format(
                config["output_location"], config["git_commit_strandphaser"]
            )
        ),
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    params:
        version=config["git_commit_strandphaser"],
        repo=config["git_repo_strandphaser"],
    shell:
        "LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {output} 2>&1"




rule config_run_summary:
    input:
        labels = "{output_folder}/config/{sample}/labels.tsv",
        info_raw="{output_folder}/counts/{sample}/{sample}.info_raw",
        ploidy_summary="{output_folder}/ploidy/{sample}/ploidy_summary.txt",
        single_paired_end_detect="{output_folder}/config/{sample}/single_paired_end_detection.txt",
    output:
        summary = "{output_folder}/config/{sample}/run_summary.txt"
    log:
        "{output_folder}/log/config/{sample}/config_run_summary.txt"
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/run_summary.py"
