#!/bin/bash

REF="/g/korbel/shared/datasets/refgenomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

snakemake \
    -j 100 \
    --configfile Snake.config-singularity.json \
    --use-singularity \
	--singularity-args "-B ${REF}:/reference.fa:ro -B ${REF}.fai:/reference.fa.fai:ro" \
	--cluster-config cluster.json \
    --local-cores 8 \
    --cluster "sbatch -o slurm/{rule}.%j.log -e slurm/{rule}.%j.log --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}" \
    --latency-wait 60 \
    --timestamp \
    --keep-going \
    --restart-times 1

