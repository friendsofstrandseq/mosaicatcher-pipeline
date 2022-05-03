rule compile_mosaic:
    """
    rule fct:
    input:
    output:
    """
    input: 
        "scripts/Mosaic/src/CMakeLists.txt"
    output:
        "scripts/Mosaic/build/mosaic"
    conda:
        "../envs/cpp_copy.yaml"
        # "mc-base"
    container:
        "docker://continuumio/miniconda3:4.11.0"
    shell:
        """
        pwd
        cd scripts/Mosaic/build
        cmake ../src
        export LD_LIBRARY_PATH=/g/korbel2/weber//lib
        make
        cd ../../..
        pwd
        """

rule install_rlib_strandphaser:
    output:
         check = touch(config['output_location'] + 'strandphaser/R_setup/strandphaser_version-{}.ok'.format(config['git_commit_strandphaser']))
    log:
        'log/strandphaser/strandphaser_install.log'
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    params:
        version = config['git_commit_strandphaser'],
        repo = config['git_repo_strandphaser']
    shell:
        'LC_MEASUREMENT=C LC_CTYPE=C TAR=$(which tar) Rscript scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo}  > {log} 2>&1'
