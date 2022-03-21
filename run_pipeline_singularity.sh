#!/bin/bash

# Set these two paths to link large external data (reference genomes) to the respective places within the image
REF="/gstock/biolo_datasets/variation/benchmark/Databases/RefSeq/ncbi-genomes-2022-02-17/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"
R_REF="/home/weber/.conda/envs/strandseqnation/lib/R/library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit"

snakemake \
    -j 6 \
    --configfile Snake.config-singularity.json \
    --use-singularity \
    --singularity-args "-B ${REF}:/reference.fa:ro \
                        -B ${REF}.fai:/reference.fa.fai:ro \
                        -B ${R_REF}:/usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit:ro" \
    --latency-wait 60 \
    --printshellcmd



