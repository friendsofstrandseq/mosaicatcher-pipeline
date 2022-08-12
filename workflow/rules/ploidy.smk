
################################################################################
# Ploidy estimation                                                            #
################################################################################

if int(config["window"]) in [50000, 100000, 200000]:



    rule estimate_ploidy:
        input:
            counts = "{output_folder}/counts/{sample}/{sample}.txt.gz"
        output:
            "{output_folder}/ploidy/{sample}/ploidy_detailled.txt"
        log:
            "{output_folder}/log/ploidy/{sample}/ploidy_detailled.log"
        threads: 48
        resources:
            mem_mb=get_mem_mb,
        params:
            # TODO move this to config
            merge_window = 1000000,
            shift_step = 500000,
            boundary_alpha = 0.05,
            max_ploidy = 6,
            add_bg_component = 'to_be_done'
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            python workflow/scripts/utils/ploidy_estimator.py --debug \
                --merge-bins-to {params.merge_window} \
                --shift-window-by {params.shift_step} \
                --max-ploidy {params.max_ploidy} \
                --boundary-alpha {params.boundary_alpha} \
                --jobs {threads} \
                --input {input.counts} \
                --output {output}
            """
            # exec = '../scripts/utils/ploidy-estimator.py --debug'
            # exec += ' --merge-bins-to {params.merge_window}'
            # exec += ' --shift-window-by {params.shift_step}'
            # exec += ' --max-ploidy {params.max_ploidy}'
            # exec += ' --boundary-alpha {params.boundary_alpha}'
            # exec += ' --jobs {threads}'
            # exec += ' --input {input.counts}'
            # exec += ' --output {output}'
            # shell(exec)
    
    checkpoint summarise_ploidy:
        input:
            ploidy = "{output_folder}/ploidy/{sample}/ploidy_detailled.txt"
        output:
            summary = "{output_folder}/ploidy/{sample}/ploidy_summary.txt"
        log:
            "{output_folder}/ploidy/{sample}/ploidy_summary.log"
        run:
            pd.read_csv(input.ploidy, sep="\t").groupby("#chrom")["ploidy_estimate"].describe().to_csv(output.summary, sep="\t")


    rule plot_ploidy:
        input:
            ploidy_detailled="{output_folder}/ploidy/{sample}/ploidy_detailled.txt"
        output:
            report(
                "{output_folder}/plots/{sample}/ploidy/{sample}.pdf",
                category="Ploidy",
            ),
        log:
            "{output_folder}/log/plot_ploidy/{sample}.log"
        conda:
            "../envs/ashleys.yaml"
        script:
            "../scripts/plotting/ploidy_plot.py"