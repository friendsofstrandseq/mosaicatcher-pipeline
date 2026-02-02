
rule fake_package:
    output:
        touch("workflow/data/ref_genomes/log/fake_package.ok"),


rule install_BSgenome_package:
    input:
        package=bsgenome_install,
    output:
        touch(
            "workflow/data/ref_genomes/config/BSgenome_{}.ok".format(
                config["reference"]
            )
        ),
    log:
        "workflow/data/ref_genomes/log/install_BSgenome_package_{}.log".format(
            config["reference"]
        ),
    params:
        selected_package=lambda wc, input: (
            "BSgenome.{}.UCSC.{}".format(
                "Mmusculus" if config["reference"] == "mm10" else "Hsapiens",
                config["reference"],
            )
            if config["reference"] in ["hg38", "hg19", "mm10"]
            else input.package
        ),
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        "Rscript workflow/scripts/utils/install_R_package.R {params.selected_package}"


rule install_sctrip_multiplot_package:
    input:
        package=bsgenome_install,
    output:
        touch("workflow/config/scTRIP_multiplot.ok"),
    log:
        "log/install_sctrip_multiplot_package.log",
    params:
        selected_package="workflow/scripts/plotting/scTRIP_multiplot/scTRIPmultiplot",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    shell:
        "Rscript workflow/scripts/utils/install_R_package.R {params.selected_package}"


rule config_run_summary:
    input:
        unpack(get_config_run_summary_inputs),
    output:
        summary=report(
            "{folder}/{sample}/config/run_summary.txt",
            category="Run summary",
            labels={"Sample": "{sample}"},
        ),
    log:
        "{folder}/log/config/{sample}/config_run_summary.txt",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/run_summary.py"


# rule install_T2T_BSgenome_tarball:
#     input:
#         tarball="workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
#     output:
#         touch("workflow/data/ref_genomes/config/BSgenome_{}.ok".format(config['reference'])),
#     params:
#         selected_package = lambda input, wc: wc.input.tarball
#     log:
#         "workflow/data/ref_genomes/log/T2T_R_tarball_install.log",
#     conda:
#         "../envs/rtools.yaml"
#     resources:
#         mem_mb=get_mem_mb_heavy,
#     script:
#         "../scripts/utils/install_R_package.R"
