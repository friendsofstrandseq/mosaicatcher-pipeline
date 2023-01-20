## Arbigent rules / Developed by Wolfram Hoeps & Hufsah Ashraf
## ---------------------------------------------------------------
## prepare_manual_segments_counts_debug: Take manual segments counts debugfile and prepare it in a way that we can easily extract the information for n informative bins per segment.
## run_regenotypeR_samplewise_singlecell: Invoke regenotype.R for each sample, creating a sv_calls_bulk and associated plots for each sample

if config["arbigent"] is True:

    rule create_hdf_file:
        input:
            mapping_track=config["arbigent_data"]["arbigent_mapability_track"],
        output:
            config["arbigent_data"]["arbigent_mapability_track_h5"],
        log:
            "workflow/data/arbigent/log/create_hdf_file.log",
        conda:
            "../envs/mc_base.yaml"
        script:
            "../scripts/arbigent_utils/create_hdf.py"

    rule custom_manual_segmentation_file:
        input:
            config["arbigent_bed_file"],
        output:
            "{folder}/{sample}/arbigent/manual_segmentation_custom.bed",
        log:
            "{folder}/{sample}/log/manual_segmentation_custom/{sample}.log",
        params:
            chromosomes="|".join(config["chromosomes"]),
        conda:
            "../envs/mc_base.yaml"
        shell:
            "grep -E -- '{params.chromosomes}' {input} > {output}"

    rule watson_crick_counts:
        input:
            bam_cells=selected_input_bam,
            bed=config["arbigent_bed_file"]
            if len(config["chromosomes"]) == 24
            else "{folder}/{sample}/arbigent/manual_segmentation_custom.bed",
            mapping=config["arbigent_data"]["arbigent_mapability_track"],
            mapping_h5=config["arbigent_data"]["arbigent_mapability_track_h5"],
        output:
            processing_counts="{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt.raw",
            debug="{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt.raw.debug",
            norm_plot_output="{folder}/{sample}/arbigent/arbigent_counts/blub.txt",
        log:
            "{folder}/log/arbigent/watson_crick_counts/{sample}/watson_crick_counts.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            bam_folder="{folder}/{sample}/selected",
            genome_chromosome_param="genome"
            if len(config["chromosomes"]) == 24
            else ",".join(config["chromosomes"]),
        threads: 16
        script:
            "../scripts/arbigent_utils/watson_crick.py"

    rule correct_watson_crick_counts:
        input:
            "{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt.raw",
        output:
            "{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt",
        log:
            "{folder}/log/arbigent/correct_watson_crick_counts/{sample}/watson_crick_counts.log",
        conda:
            "../envs/mc_base.yaml"
        shell:
            "sed 's/.sort.mdup//g' {input} > {output}"

    rule prepare_manual_segments_counts_debug:
        input:
            counts_file="{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt.raw.debug",
        output:
            msc="{folder}/{sample}/arbigent/sv_calls/msc.debug",
        log:
            "{folder}/log/arbigent/sv_calls/{sample}/watson_crick_counts.log",
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            awk '!seen[$1,$2,$3]++' {input.counts_file} > {output.msc}
            """

    rule mosaiClassifier_calc_probs_arbigent:
        input:
            # [W] windows_specs is x_fixed_norm. 
            counts="{folder}/{sample}/counts/{sample}.txt.gz",
            info="{folder}/{sample}/counts/{sample}.info",
            states="{folder}/{sample}/strandphaser/StrandPhaseR_final_output.txt",
            bp="{folder}/{sample}/arbigent/arbigent_counts/manual_segments_counts.txt",
        output:
            output="{folder}/{sample}/arbigent_mosaiclassifier/sv_probabilities/probabilities.Rdata",
        log:
            "{folder}/log/arbigent/mosaiClassifier_calc_probs_arbigent/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier.snakemake.R"
