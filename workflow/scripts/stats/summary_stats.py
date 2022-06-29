#!/usr/bin/python

import subprocess
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
    p.append(snakemake.input.merged)
additional_params = ' '.join(p)
subprocess.call('workflow/scripts/stats/callset_summary_stats.py --segmentation {} --strandstates {} --complex-regions {} {} {}  > {} '.format(
    snakemake.input.segmentation,
    snakemake.input.strandstates,
    snakemake.input.complex,
    additional_params,
    snakemake.input.sv_calls,
    snakemake.output.tsv
    ), shell=True
)
