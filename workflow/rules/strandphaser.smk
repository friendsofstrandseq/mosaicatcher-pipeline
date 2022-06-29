# from workflow.scripts.utils.utils import get_mem_mb 

################################################################################
# StrandPhaseR things                                                          #
################################################################################

rule convert_strandphaser_input:
    """
    rule fct: extract only segmentation with WC orientation 
    input: initial_strand_state file coming from rule segmentation_selection & info file from mosaic count output
    output: filtered TSV file with start/end coordinates of WC-orientated segment to be used by strandphaser
    """
    input:
        states = "{output}/segmentation/{sample}/Selection_initial_strand_state",
        info   = "{output}/counts/{sample}/{sample}.info"
    output:
        "{output}/strandphaser/{sample}/strandphaser_input.txt"
    log:
        "{output}/log/strandphaser/convert_strandphaser_input/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    script:
        "../scripts/strandphaser_scripts/helper.convert_strandphaser_input.R"

checkpoint determine_sex_per_cell:
    """
    rule fct:
    input:
    output:
    """
    input:
        expand("{input_folder}/{sample}/selected/", input_folder=config["input_bam_location"], sample=samples)
    output:
        sex_analysis_cellwise = "{output}/config/{sample}/sex_analysis_cells.tsv",
        sex_analysis_samplewise = "{output}/config/{sample}/sex_analysis_sample.txt"
    log:
        "{output}/log/strandphaser/determine_sex_per_cell/{sample}.log"
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/chrxy_analysis.py"



# TODO : replace by clean config file if possible or by temporary removed file 
rule prepare_strandphaser_config_per_chrom:
    """
    rule fct: prepare config file used by strandphaser
    input: input used only for wildcards : sample, window & bpdens
    output: config file used by strandphaser
    """
    input:
        "{output}/segmentation/{sample}/Selection_initial_strand_state"
    output:
        "{output}/strandphaser/{sample}/StrandPhaseR.{chrom}.config"
    log:
        "{output}/log/strandphaser/{sample}/StrandPhaseR.{chrom}.log"
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/strandphaser_scripts/prepare_strandphaser.py"




rule run_strandphaser_per_chrom:
    """
    rule fct: run strandphaser for each chromosome 
    input: strandphaser_input.txt from rule convert_strandphaser_input ; genotyped snv for each chrom by freebayes ; configfile created by rule prepare_strandphaser_config_per_chrom ; bam folder
    output:
    """
    input:
        install_strandphaser = rules.install_rlib_strandphaser.output,
        wcregions    = "{output}/strandphaser/{sample}/strandphaser_input.txt",
        snppositions = "{output}/snv_genotyping/{sample}/{chrom}.vcf",
        configfile   = "{output}/strandphaser/{sample}/StrandPhaseR.{chrom}.config",
        # bamfolder    = config["input_bam_location"] + "{sample}/selected",
        # TODO : tmp solution
    output:
        "{output}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        "{output}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf",
        # "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",
        report(
            "{output}/strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",
            category="StrandPhaseR",
            subcategory = "{sample}",
            caption="../report/strandphaser_haplotypes.rst",
            labels={"Sample" : "{sample}", "Chrom" : "{chrom}"}
        )
    log:
        "{output}/log/run_strandphaser_per_chrom/{sample}/{chrom}.log"
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb = get_mem_mb,
    params:
        input_bam = lambda wc: "{}/{}/selected".format(config["input_bam_location"], wc.sample),
        output = lambda wc: "{}/strandphaser/{}/StrandPhaseR_analysis.{}".format(config["output_location"], wc.sample, wc.chrom)
    shell:
        # {config[Rscript]}
        """
        Rscript workflow/scripts/strandphaser_scripts/StrandPhaseR_pipeline.R \
                {params.input_bam} \
                {params.output} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ > {log} 2>&1
        """

rule merge_strandphaser_vcfs:
    input:
        ## OLD calling
        # vcfs=expand("strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"]),
        # tbis=expand("strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"]),
        ## NEW calling that takes into account sex sample (see checkpoint determine_sex_per_cell)
        vcfs= ancient(aggregate_vcf_gz),
        tbis= ancient(aggregate_vcf_gz_tbi)
    output:
        vcfgz="{output}/strandphaser/phased-snvs/{sample}.vcf.gz"
    log:
        "{output}/log/merge_strandphaser_vcfs/{sample}.log"
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb = get_mem_mb,
    shell:
        """
        (bcftools concat -a {input.vcfs} | bcftools view -o {output.vcfgz} -O z --genotype het --types snps - ) > {log} 2>&1
        """



rule combine_strandphaser_output:
    input:
        aggregate_phased_haps
    output:
         "{output}/strandphaser/{sample}/strandphaser_phased_haps_merged.txt"
    log:
        "{output}/log/combine_strandphaser_output/{sample}.log"
    resources:
        mem_mb = get_mem_mb,
    run:
        ## Errors on slurm with previous version
        # """
        # set -o pipefail;
        # cat {input} | head -n1 > {output};
        # tail -q -n+2 {input} >> {output}; > {log} 2>&1
        # """
        ## New version using pandas
        import pandas as pd
        pd.concat([pd.read_csv(file, sep="\t") for j, file in enumerate(input)]).to_csv(output[0], sep="\t", index=False)
        


rule convert_strandphaser_output:
    input:
        phased_states  = "{output}/strandphaser/{sample}/strandphaser_phased_haps_merged.txt",
        initial_states = "{output}/segmentation/{sample}/Selection_initial_strand_state",
        info           = "{output}/counts/{sample}/{sample}.info"
    output:
        "{output}/strandphaser/{sample}/StrandPhaseR_final_output.txt"
    log:
        "{output}/log/convert_strandphaser_output/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb = get_mem_mb,
    script:
        "../scripts/strandphaser_scripts/helper.convert_strandphaser_output.R"
