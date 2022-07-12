from workflow.scripts.utils.utils import get_mem_mb 

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
        states = config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state",
        info   = config["output_location"] + "counts/{sample}/{sample}.info"
    output:
        config["output_location"] + "strandphaser/{sample}/strandphaser_input.txt"
    log:
        config["output_location"] + "log/strandphaser/convert_strandphaser_input/{sample}.log"
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
        config["input_bam_location"] + "{sample}/selected/"
    output:
        sex_analysis_cellwise = config["output_location"] + "config/{sample}/sex_analysis_cells.tsv",
        sex_analysis_samplewise = config["output_location"] + "config/{sample}/sex_analysis_sample.txt"
    log:
        config["output_location"] + "log/strandphaser/determine_sex_per_cell/{sample}.log"
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/utils/chrxy_analysis.py"


def aggregate_phased_haps(wildcards):
    with checkpoints.determine_sex_per_cell.get(sample=wildcards.sample).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split('\t')[1]
        if sex == "M":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrX", "chrY"]]
        elif sex == "F":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrY"]]
        return expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt", chrom=config["chromosomes"])



def aggregate_vcf_gz(wildcards):
    with checkpoints.determine_sex_per_cell.get(sample=wildcards.sample).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split('\t')[1]
        if sex == "M":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrX", "chrY"]]
        elif sex == "F":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrY"]]
        return expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"])



def aggregate_vcf_gz_tbi(wildcards):
    with checkpoints.determine_sex_per_cell.get(sample=wildcards.sample).output.sex_analysis_samplewise.open() as f:
        sex = f.read().strip().split('\t')[1]
        if sex == "M":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrX", "chrY"]]
        elif sex == "F":
            config["chromosomes"] = [c for c in config["chromosomes"] if c not in ["chrY"]]
        return expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"])



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
        # config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",
        report(
            config["output_location"] + "strandphaser/{sample}/StrandPhaseR_analysis.{chrom}/SingleCellHaps/{chrom}_singleCellHaps.pdf",
            category="StrandPhaseR",
            subcategory = "{sample}",
            caption="../report/strandphaser_haplotypes.rst",
            labels={"Sample" : "{sample}", "Chrom" : "{chrom}"}
        )
    log:
        config["output_location"] + "log/run_strandphaser_per_chrom/{sample}/{chrom}.log"
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb = get_mem_mb,
        time = "02:00:00",
    shell:
        # {config[Rscript]}
        """
        Rscript workflow/scripts/strandphaser_scripts/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                {config[output_location]}strandphaser/{wildcards.sample}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ > {log} 2>&1
        """

rule merge_strandphaser_vcfs:
    input:
        ## OLD calling
        # vcfs=expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz", chrom=config["chromosomes"]),
        # tbis=expand(config["output_location"] + "strandphaser/{{sample}}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf.gz.tbi", chrom=config["chromosomes"]),
        ## NEW calling that takes into account sex sample (see checkpoint determine_sex_per_cell)
        vcfs= ancient(aggregate_vcf_gz),
        tbis= ancient(aggregate_vcf_gz_tbi)
    output:
        vcfgz=config["output_location"] + "strandphaser/phased-snvs/{sample}.vcf.gz"
    log:
        config["output_location"] + "log/merge_strandphaser_vcfs/{sample}.log"
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
        config["output_location"] +  "strandphaser/{sample}/strandphaser_phased_haps_merged.txt"
    log:
        config["output_location"] + "log/combine_strandphaser_output/{sample}.log"
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
        phased_states  = config["output_location"] + "strandphaser/{sample}/strandphaser_phased_haps_merged.txt",
        initial_states = config["output_location"] + "segmentation/{sample}/Selection_initial_strand_state",
        info           = config["output_location"] + "counts/{sample}/{sample}.info"
    output:
        config["output_location"] + "strandphaser/{sample}/StrandPhaseR_final_output.txt"
    log:
        config["output_location"] + "log/convert_strandphaser_output/{sample}.log"
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb = get_mem_mb,
    script:
        "../scripts/strandphaser_scripts/helper.convert_strandphaser_output.R"
