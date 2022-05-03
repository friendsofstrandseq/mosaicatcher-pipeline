
################################################################################
# StrandPhaseR things                                                          #
################################################################################


# TODO : replace R script by integrating directly pandas in the pipeline / potentialy use piped output to following rule ?

rule convert_strandphaser_input:
    """
    rule fct: extract only segmentation with WC orientation 
    input: initial_strand_state file coming from rule segmentation_selection & info file from mosaic count output
    output: filtered TSV file with start/end coordinates of WC-orientated segment to be used by strandphaser
    """
    input:
        states = config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state",
        # URGENT : hard coded 500000 file name ???
        # info   = config["output_location"] + "counts/{sample}/500000.info"
        # FIXME : quick workaround with {window} wc
        info   = config["output_location"] + "counts/{sample}/{sample}.info"
    output:
        config["output_location"] + "strandphaser/{sample}/strandphaser_input.txt"
    log:
        config["output_location"] + "log/strandphaser/convert_strandphaser_input/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/strandphaser_scripts/helper.convert_strandphaser_input.R"


# TODO : make something similar to mosaic with C++ dep
# CHECKME : check if possible to write something more snakemak"ic" & compliant with conda/singularity running env
# WARNING : I/O path definition
# WARNING : Try to find a solution to install stranphaser in a conda environment => contact david porubsky to move on the bioconductor ?

# rule install_StrandPhaseR:
#     output:
#         "utils/R-packages/StrandPhaseR/R/StrandPhaseR"
#     log:
#         "log/install_StrandPhaseR.log"
#     conda:
#         "../envs/strandphaser.yaml"
#     shell:
#         """
#         Rscript utils/install_strandphaser.R > {log} 2>&1
#         """

# ruleorder: install_rlib_strandphaser > run_strandphaser_per_chrom

# rule install_rlib_strandphaser:
#     output:
#          check = touch(config['output_location'] + 'strandphaser/R_setup/strandphaser_version-{}.ok'.format(config['git_commit_strandphaser']))
#     log:
#         config["output_location"] + 'log/strandphaser/strandphaser_install.log'
#     conda:
#         "../envs/rtools.yaml"
#     resources:
#         mem_total_mb = 4096,
#         mem_per_cpu_mb = 4096
#     params:
#         version = config['git_commit_strandphaser'],
#         repo = config['git_repo_strandphaser']
#     shell:
#         'LC_MEASUREMENT=C LC_CTYPE=C TAR=$(which tar) Rscript scripts/strandphaser_scripts/install_strandphaser.R {params.version} {params.repo} &> {log}'


# TODO : replace by clean config file if possible or by temporary removed file 
rule prepare_strandphaser_config_per_chrom:
    """
    rule fct: prepare config file used by strandphaser
    input: input used only for wildcards : sample, window & bpdens
    output: config file used by strandphaser
    """
    input:
        config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state"
    output:
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR.{chrom}.config"
    run:
        with open(output[0], "w") as f:
            print("[General]",                    file = f)
            print("numCPU           = 1",         file = f)
            print("chromosomes      = '" + wildcards.chrom + "'", file = f)
            if (config["paired_end"]):
                print("pairedEndReads   = TRUE",  file = f)
            else:
                print("pairedEndReads   = FALSE", file = f)
            print("min.mapq         = 10",        file = f)
            print("",                             file = f)
            print("[StrandPhaseR]",               file = f)
            print("positions        = NULL",      file = f)
            print("WCregions        = NULL",      file = f)
            print("min.baseq        = 20",       file = f)
            print("num.iterations   = 2",        file = f)
            print("translateBases   = TRUE",     file = f)
            print("fillMissAllele   = NULL",     file = f)
            print("splitPhasedReads = TRUE",     file = f)
            print("compareSingleCells = TRUE",     file = f)
            print("callBreaks       = FALSE",    file = f)
            print("exportVCF        = '", wildcards.sample, "'", sep = "", file = f)
            print("bsGenome         = '", config["R_reference"], "'", sep = "", file = f)




rule run_strandphaser_per_chrom:
    """
    rule fct: run strandphaser for each chromosome 
    input: strandphaser_input.txt from rule convert_strandphaser_input ; genotyped snv for each chrom by freebayes ; configfile created by rule prepare_strandphaser_config_per_chrom ; bam folder
    output:
    """
    input:
        wcregions    = config["output_location"] + "strandphaser/{sample}/strandphaser_input.txt",
        snppositions = config["output_location"] + "snv_genotyping/{sample}/{chrom}.vcf",
        configfile   = config["output_location"] + "strandphaser/{sample}/StrandPhaseR.{chrom}.config",
        bamfolder    = config["input_bam_location"] + "{sample}/selected",
        # TODO : tmp solution
        strandphaser_install = config['output_location'] + 'strandphaser/R_setup/strandphaser_version-{}.ok'.format(config['git_commit_strandphaser'])
    output:
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf",
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",

        # report(
        #     config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",
        #     category="StrandPhaseR",
        #     subcategory = "{sample}",
        #     caption="../report/strandphaser_haplotypes.rst",
        #     labels={"Sample" : "{sample}", "Chrom" : "{chrom}"}
        # )
    log:
        "log/run_strandphaser_per_chrom/{sample}/{chrom}.log"
    conda:
        "../envs/rtools.yaml"
    shell:
        # {config[Rscript]}
        """
        Rscript scripts/strandphaser_scripts/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                {config[output_location]}strandphaser/{wildcards.sample}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \

        """

rule merge_strandphaser_vcfs:
    input:
        vcfs=expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"]),
        tbis=expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"]),
    output:
        vcfgz=config["output_location"] + "strandphaser/phased-snvs/{sample}.vcf.gz"
    log:
        "log/merge_strandphaser_vcfs/{sample}.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        """
        (bcftools concat -a {input.vcfs} | bcftools view -o {output.vcfgz} -O z --genotype het --types snps - ) > {log} 2>&1
        """


rule combine_strandphaser_output:
    input:
        expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt", chrom = config["chromosomes"])
    output:
        config["output_location"] +  "strandphaser/{sample}/strandphaser_phased_haps_merged.txt"
    log:
        "log/combine_strandphaser_output/{sample}.log"
    shell:
        """
        set +o pipefail
        cat {input} | head -n1 > {output};
        tail -q -n+2 {input} >> {output};
        """


rule convert_strandphaser_output:
    input:
        phased_states  = config["output_location"] + "strandphaser/{sample}/strandphaser_phased_haps_merged.txt",
        initial_states = config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state",
        # info           = config["output_location"] + "counts/{sample}/500000_fixed.info"
        info           = config["output_location"] + "counts/{sample}/{sample}.info"
    output:
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR_final_output.txt"
    log:
        "log/convert_strandphaser_output/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/strandphaser_scripts/helper.convert_strandphaser_output.R"

