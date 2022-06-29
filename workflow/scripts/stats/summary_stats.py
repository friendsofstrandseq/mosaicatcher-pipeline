#!/usr/bin/python

p = []
try:
    f = snakemake.config["ground_truth_clonal"][snakemake.wildcards.sample]
    if len(f) > 0:
        p.append('--true-events-clonal')
        p.append(f)
except KeyError:
    pass
try:
    f = snakemake.config["ground_truth_single_cell"][snakemake.wildcards.sample]
    if len(f) > 0:
        p.append('--true-events-single-cell')
        p.append(f)
except KeyError:
    pass
if snakemake.wildcards.filter == 'TRUE':
    p.append('--merged-file')
    p.append(input.merged)
additional_params = ' '.join(p)
shell('workflow/scripts/stats/callset_summary_stats.py --segmentation {snakemake.input.segmentation} --strandstates {snakemake.input.strandstates} --complex-regions {input.complex} {additional_params} {snakemake.input.sv_calls}  > {snakemake.output.tsv} ')
