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
        mem_mb="10G",
    params:
        version=config["git_commit_strandphaser"],
        repo=config["git_repo_strandphaser"],
    shell:
        "LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {output} 2>&1"
