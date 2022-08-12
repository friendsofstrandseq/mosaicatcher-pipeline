# rule install_ashleys:
#     output:
#         touch("{output_folder}/ashleys/ashleys_install_success.txt"),
#     log:
#         "{output_folder}/log/ashleys/install_ashleys.log",
#     conda:
#         "../envs/ashleys.yaml"
#     shell:
#         "( git clone https://github.com/friendsofstrandseq/ashleys-qc.git &&"
#         "cd ashleys-qc &&"
#         "python setup.py install ) &> {log}"


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
        mem_mb=8000,
    params:
        version=config["git_commit_strandphaser"],
        repo=config["git_repo_strandphaser"],
    shell:
        "LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {output} 2>&1"
