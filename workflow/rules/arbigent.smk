if config["arbigent"] is True:

    ################ RULES FOR PART 1:
    ############### RUNNING REGENOTYPER MAIN PART, up to Phasing

    rule run_regenotypeR_samplewise_bulk:
        """
        Invoke regenotype.R for each sample, creating a sv_calls_bulk and associated
        plots for each sample
        """
        input:
            probabilities_table="{folder}/{sample}/arbigent/mosaiclassifier/sv_probabilities/probabilities.Rdata",
            msc="{folder}/{sample}/arbigent/sv_calls/msc.debug",
        output:
            sv_calls_bulk_dir=directory(
                "{folder}/{sample}/arbigent/regenotyper_samplewise_bulk/sv_calls/"
            ),
            sv_calls_bulk="{folder}/{sample}/arbigent/regenotyper_samplewise_bulk/sv_calls/all/sv_calls_bulk.txt",
        log:
            "{folder}/log/run_regenotypeR_samplewise_bulk/{sample}.log",
        # params:
        #     outputfolder=lambda wc, output: f"{'/'.join(str(output.sv_calls_bulk).split('/')[:-1])}/",
        conda:
            "../envs/dev/rtools_enhanced.yaml"
        shell:
            """
            Rscript workflow/scripts/arbigent/regenotype.R \
                            -f {input.probabilities_table} \
                            -c {input.msc} \
                            -o {output.sv_calls_bulk_dir}/ > {log} 2>&1
            """

    rule merge_alls:
        """
        Take the individual sv_calls_bulk created for each sample, and merge them together.
        In merging we have to take into account to only keep the header of the first file but not the others.
        This is done with awk.
        """
        input:
            sv_calls_bulk="{folder}/{sample}/arbigent/regenotyper_samplewise_bulk/sv_calls/all/sv_calls_bulk.txt",
        output:
            "{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/all_sv_calls_unphased.txt",
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            awk 'FNR==1 && NR!=1 {{ while (/^chrom/) getline; }} 1 {{print}}' {input.sv_calls_bulk} > {output}
            """

    rule rephase_all_txt:
        """
        If no phasing should be perfomed, we just copy over the old file.
        """
        input:
            all_txt="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/all_sv_calls_unphased.txt",
        output:
            all_txt_rephased="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/all_sv_calls_phased.txt",
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            cp {input.all_txt} {output.all_txt_rephased}
            """

    rule make_output_vcfs:
        """
        Take the phased sv tables, and turn them into a series of viable vcfs.
        The expand({folder}) part seems a bit unnecessary, I'm sure this can be done cleaner.
        If there is time I will revisit. But anyway it works now.
        """
        input:
            alltxt="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/all_sv_calls_phased.txt",
            # WARNING
            msc="{folder}/{sample}/arbigent/sv_calls/msc.debug",
        output:
            vcf1="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/res_all.vcf",
            res_csv="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/res.csv",
            res_csv_dir=directory(
                "{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/"
            ),
        log:
            "{folder}/log/make_output_vcfs/{sample}.log",
        params:
            outdir=lambda wc, output: f"{'/'.join(str(output.res_csv).split('/')[:-1])}/",
        conda:
            "../envs/dev/rtools_enhanced.yaml"
        shell:
            """
            Rscript workflow/scripts/arbigent/table_to_vcfs.R \
                                    -a {input.alltxt} \
                                    -m {input.msc} \
                                    -o {output.res_csv_dir}/ > {log} 2>&1
            """

    rule add_verdict:
        """
        Take res_csv, and return a similar file - res_verdicted.vcf, which contains (surprise) our verdict
        """
        input:
            res_csv="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/res.csv",
        output:
            verdicted_table="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/res_verdicted.vcf",
        log:
            "{folder}/log/add_verdict/{sample}.log",
        params:
            names_gm_to_na=1,
        conda:
            "../envs/dev/rtools_enhanced.yaml"
        shell:
            """
            Rscript workflow/scripts/arbigent/add_filter.R \
                                -i {input.res_csv} \
                                -n {params.names_gm_to_na} \
                                -o {output.verdicted_table} > {log} 2>&1
            """

    rule qc_result:
        input:
            verdicted_table="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/arbigent_results/res_verdicted.vcf",
        output:
            verdict_plot="{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/qc/lineplot_gts.pdf",
            verdict_plot_dir=directory(
                "{folder}/{sample}/arbigent/regenotyper_allsamples_bulk/qc/"
            ),
        log:
            "{folder}/log/qc_result/{sample}.log",
        params:
            qc_folder=lambda wc, output: f"{'/'.join(str(output.verdict_plot).split('/')[:-1])}/",
        conda:
            "../envs/dev/rtools_enhanced.yaml"
        shell:
            """
            Rscript workflow/scripts/arbigent/qc_res_verdicted.R \
                                -f {input.verdicted_table} \
                                -o {output.verdict_plot_dir}/ > {log} 2>&1
            """
