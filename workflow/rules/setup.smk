rule install_rlib_strandphaser:
    output:
         check = touch(config['output_location'] + 'strandphaser/R_setup/strandphaser_version-{}.ok'.format(config['git_commit_strandphaser']))
    log:
        config['output_location'] + 'log/strandphaser/strandphaser_install.log'
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb = "10G",
    params:
        version = config['git_commit_strandphaser'],
        repo = config['git_repo_strandphaser']
    shell:
        'LC_CTYPE=C TAR=$(which tar) Rscript workflow/scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {log} 2>&1'
