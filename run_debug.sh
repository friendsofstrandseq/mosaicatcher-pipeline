#!/bin/bash

~/.pixi/bin/pixi run snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --sdm conda --conda-frontend mamba --nolock \
    --force .tests/data_CHR17/RPE-BM510/plots/sv_clustering/stringent-filterTRUE-chromosome.pdf \
    --skip-script-cleanup -p
