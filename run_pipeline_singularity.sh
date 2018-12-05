#!/bin/bash

# Set these two paths to link large external data (reference genomes) to the respective places within the image
REF="/g/korbel/shared/datasets/refgenomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
R_REF="/g/korbel/meiers/R-lib/3.4.0-foss-2016b/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit"

snakemake \
    -j 6 \
    --configfile Snake.config-singularity.json \
    --use-singularity \
    --singularity-args "-B ${REF}:/reference.fa:ro \
                        -B ${REF}.fai:/reference.fa.fai:ro \
                        -B ${R_REF}:/usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit:ro" \
    --latency-wait 60 \
    --printshellcmd

