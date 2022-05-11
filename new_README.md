Prerequire : singularity + conda 

conda install -c conda-forge mamba
mamba create -n mosaicatcher_env -c bioconda -c conda-forge snakemake pandas

git clone https://git.embl.de/tweber/mosaicatcher-update.git

snakemake \
    --use-conda  \
    --cores 40  \
    --config plot=True \
    mode=mosaiclassifier \
    output_location=/g/korbel2/weber/Mosaicatcher_output/Mosaicatcher_output_singularity_test_TALL \
    input_bam_location=/g/korbel2/weber/MosaiCatcher_files/T-ALL_P1  \
    -p \
    --conda-frontend mamba \
    --use-singularity \
    --singularity-args "-B /g:/g"