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
<<<<<<< HEAD
        check=touch(
            expand(
                "{}/strandphaser/R_setup/strandphaser_version-{}.ok".format(
                    config["output_location"], config["git_commit_strandphaser"]
                )
            )
        ),
=======
        check = touch(config['output_location'] + 'strandphaser/R_setup/strandphaser_version-{}.ok'.format(config['git_commit_strandphaser']))
>>>>>>> master
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
<<<<<<< HEAD
        "LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {output} 2>&1"
=======
        'LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {log} 2>&1'

rule cp_yaml_config:
    input:
        "config/config.yaml"
    output:
        config["output_location"] + "config/config.yaml"
    shell:
        "cp {input} {output}"
>>>>>>> master
