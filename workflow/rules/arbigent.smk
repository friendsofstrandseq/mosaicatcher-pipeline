if config["arbigent"] is True:

    ################ RULES FOR PART 1:
    ############### RUNNING REGENOTYPER MAIN PART, up to Phasing

    rule merge_alls:
"""
        Take the individual sv_calls_bulk created for each sample, and merge them together. 
        In merging we have to take into account to only keep the header of the first file but not the others.
        This is done with awk.
        """
        input:
            sv_calls_bulk=expand(
                "{path}/regenotyper_samplewise_{{mode}}/{sample}/all/sv_calls_bulk.txt",
                path=path_to_pipe,
                sample=SAMPLES,
            ),
        output:
            expand(
                "{path}/regenotyper_allsamples_{{mode}}/all_sv_calls_unphased.txt",
                path=path_to_pipe,
            ),
        shell:
            """
            awk 'FNR==1 && NR!=1 {{ while (/^chrom/) getline; }} 1 {{print}}' {input.sv_calls_bulk} > {output}
            """

    rule rephase_all_txt:
"""
        If no phasing should be perfomed, we just copy over the old file.
        """
        input:
            all_txt=expand(
                "{path}/regenotyper_allsamples_{{mode}}/all_sv_calls_unphased.txt",
                path=path_to_pipe,
            ),
        output:
            all_txt_rephased=expand(
                "{path}/regenotyper_allsamples_{{mode}}/all_sv_calls_phased.txt",
                path=path_to_pipe,
            ),
        shell:
            """
            cp {input.all_txt} {output.all_txt_rephased}
            """

    rule make_output_vcfs:
"""
        Take the phased sv tables, and turn them into a series of viable vcfs.
        The expand({path}) part seems a bit unnecessary, I'm sure this can be done cleaner. 
        If there is time I will revisit. But anyway it works now. 
        """
        input:
            alltxt=expand(
                "{path}/regenotyper_allsamples_{{mode}}/all_sv_calls_phased.txt",
                path=path_to_pipe,
            ),
            msc=expand(
                "{path}/regenotyper_samplewise_{{mode}}/msc.debug", path=path_to_pipe
            ),
        output:
            vcf1=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/res_all.vcf",
                path=path_to_pipe,
            ),
            res_csv=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/res.csv",
                path=path_to_pipe,
            ),
        params:
            outdir=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results",
                path=path_to_pipe,
            ),
        shell:
            """
            {rscript_command} table_to_vcfs.R \
                                    -a {input.alltxt} \
                                    -m {input.msc} \
                                    -o {params.outdir}
            """

    rule add_verdict:
"""
        Take res_csv, and return a similar file - res_verdicted.vcf, which contains (surprise) our verdict
        """
        input:
            res_csv=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/res.csv",
                path=path_to_pipe,
            ),
        output:
            verdicted_table=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/res_verdicted.vcf",
                path=path_to_pipe,
            ),
        params:
            names_gm_to_na=1,
        shell:
            """
            {rscript_command} add_filter.R \
                                -i {input.res_csv} \
                                -n {params.names_gm_to_na} \
                                -o {output.verdicted_table}
            """

    rule qc_result:
        input:
            verdicted_table=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/res_verdicted.vcf",
                path=path_to_pipe,
            ),
        output:
            verdict_plot=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/qc/lineplot_gts.pdf",
                path=path_to_pipe,
            ),
        params:
            qc_folder=expand(
                "{path}/regenotyper_allsamples_{{mode}}/arbigent_results/qc",
                path=path_to_pipe,
            ),
        shell:
            """
            {rscript_command} qc_res_verdicted.R \
                                -f {input.verdicted_table} \
                                -o {params.qc_folder}
            """
