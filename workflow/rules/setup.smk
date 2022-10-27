# rule install_rlib_strandphaser:
#     output:
#         check=touch(
#             expand(
#                 "{folder}/config/strandphaser/R_setup/strandphaser_version-{commit}.ok".format(
#                     folder=config["data_location"],
#                     commit=config["git_commit_strandphaser"],
#                 )
#             )
#         ),
#     log:
#         expand(
#             "{}/log/strandphaser/R_setup/strandphaser_version-{}.ok".format(
#                 config["data_location"], config["git_commit_strandphaser"]
#             )
#         ),
#     conda:
#         "../envs/rtools.yaml"
#     resources:
#         mem_mb=get_mem_mb_heavy,
#     params:
#         version=config["git_commit_strandphaser"],
#         repo=config["git_repo_strandphaser"],
#     shell:
#         "LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {output} 2>&1"


rule config_run_summary:
    input:
        labels="{folder}/{sample}/config/labels.tsv",
        info_raw="{folder}/{sample}/counts/{sample}.info_raw",
        ploidy_summary="{folder}/{sample}/ploidy/ploidy_summary.txt",
        single_paired_end_detect="{folder}/{sample}/config/single_paired_end_detection.txt",
    output:
        summary=report("{folder}/config/{sample}/run_summary.txt", category="Run summary", labels={"Sample": "{sample}"}),
    log:
        "{folder}/log/config/{sample}/config_run_summary.txt",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/run_summary.py"
