import os

# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

# HTTP = HTTPRemoteProvider()


rule dl_example_data:
    input:
        storage(
            "https://sandbox.zenodo.org/record/1074721/files/TEST_EXAMPLE_DATA.zip",
            keep_local=True,
        ),
    output:
        touch("config/dl_example_data.ok"),
    log:
        touch("log/config/dl_example_data.ok"),
    run:
        shell("unzip {input} -d .")


rule download_hg19_reference:
    input:
        storage(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg19.fa",
    log:
        "workflow/data/ref_genomes/log/hg19.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg19.fa.gz
        gunzip workflow/data/ref_genomes/hg19.fa.gz
        """


rule download_hg38_reference:
    input:
        storage(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg38.fa",
    log:
        "workflow/data/ref_genomes/log/hg38.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg38.fa.gz
        gunzip workflow/data/ref_genomes/hg38.fa.gz
        """


rule download_T2T_reference:
    input:
        storage(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/T2T.fa",
    log:
        "workflow/data/ref_genomes/log/T2T.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/T2T.fa.gz
        gunzip workflow/data/ref_genomes/T2T.fa.gz
        """


rule download_mm10_reference:
    input:
        storage(
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/mm10.fa",
    log:
        "workflow/data/ref_genomes/log/mm10.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/mm10.fa.gz
        gunzip workflow/data/ref_genomes/mm10.fa.gz
        """


rule download_T2T_tarball:
    input:
        storage(
            "https://zenodo.org/record/7697400/files/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
    log:
        "workflow/data/ref_genomes/log/T2T_tarball.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz
        """


rule download_arbigent_mappability_track:
    input:
        storage(
            "https://zenodo.org/record/7697400/files/mapping_counts_allchrs_hg38.txt",
            keep_local=True,
        ),
    output:
        config["arbigent_data"]["arbigent_mapability_track"],
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    shell:
        """
        directory="workflow/data/arbigent/"
        mkdir -p "$directory"
        mv {input} {output}
        """


rule download_scnova_data:
    input:
        ancient(
            storage(
                "https://zenodo.org/record/7697400/files/scNOVA_data_models.zip",
                keep_local=True,
            )
        ),
    output:
        "workflow/data/scNOVA/utils/bin_chr_length.bed",
        "workflow/data/scNOVA/utils/bin_Genebody_all.bed",
        "workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort_ann_sort_GC_ensemble.txt",
        "workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort.txt",
        "workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
        "workflow/data/scNOVA/utils/bin_Genes_for_CNN_sort.txt.corrected",
        "workflow/data/scNOVA/utils/Deeptool_Genes_for_CNN_merge_sort_lab_final.txt",
        "workflow/data/scNOVA/utils/Features_reshape_CpG_orientation_impute.txt",
        "workflow/data/scNOVA/utils/Features_reshape_CpG_orientation.txt",
        "workflow/data/scNOVA/utils/Features_reshape_GC_orientation_impute.txt",
        "workflow/data/scNOVA/utils/Features_reshape_GC_orientation.txt",
        "workflow/data/scNOVA/utils/Features_reshape_RT_orientation.txt",
        "workflow/data/scNOVA/utils/Features_reshape_size_orientation.txt",
        "workflow/data/scNOVA/utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
        "workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
        "workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort.bed",
        "workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
        "workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SVM.txt",
        "workflow/data/scNOVA/utils/Strand_seq_matrix_TES_for_SVM.txt",
        "workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    conda:
        "../envs/scNOVA/scNOVA_DL.yaml"
    # container:
    #     None
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/scNOVA_data_models.zip
        unzip workflow/data/scNOVA_data_models.zip -d workflow/data/
        """
