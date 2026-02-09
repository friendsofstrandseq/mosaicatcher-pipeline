# Register HTTP storage provider for downloading reference genomes
storage http:
    provider="http",
    max_requests_per_second=10,


# Get reference base directory from config (supports multi-user HPC setups)
REF_BASE_DIR = config.get("reference_base_dir", "workflow/data/ref_genomes")


rule dl_example_data:
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
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
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/hg19.fa",
    log:
        f"{REF_BASE_DIR}/log/hg19.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/hg19.fa.gz
        gunzip {REF_BASE_DIR}/hg19.fa.gz
        """


rule download_hg38_reference:
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/hg38.fa",
    log:
        f"{REF_BASE_DIR}/log/hg38.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/hg38.fa.gz
        gunzip {REF_BASE_DIR}/hg38.fa.gz
        """


rule download_T2T_reference:
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/T2T.fa",
    log:
        f"{REF_BASE_DIR}/log/T2T.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/T2T.fa.gz
        gunzip {REF_BASE_DIR}/T2T.fa.gz
        """


rule download_mm10_reference:
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/mm10.fa",
    log:
        f"{REF_BASE_DIR}/log/mm10.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/mm10.fa.gz
        gunzip {REF_BASE_DIR}/mm10.fa.gz
        """


rule download_mm39_reference:
    localrule: True
    input:
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/mm39.fa",
    log:
        f"{REF_BASE_DIR}/log/mm39.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/mm39.fa.gz
        gunzip {REF_BASE_DIR}/mm39.fa.gz
        """


rule download_canFam3_reference:
    localrule: True
    input:
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/canFam3.fa",
    log:
        f"{REF_BASE_DIR}/log/canFam3.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/canFam3.fa.gz
        gunzip {REF_BASE_DIR}/canFam3.fa.gz
        """


rule download_canFam4_reference:
    localrule: True
    input:
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.fa.gz",
            keep_local=True,
        ),
    output:
        f"{REF_BASE_DIR}/canFam4.fa",
    log:
        f"{REF_BASE_DIR}/log/canFam4.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        f"""
        directory="{REF_BASE_DIR}/"
        mkdir -p "$directory" "$directory/log"
        mv {{input}} {REF_BASE_DIR}/canFam4.fa.gz
        gunzip {REF_BASE_DIR}/canFam4.fa.gz
        """


rule generate_canFam3_bin_bed:
    input:
        fasta="workflow/data/ref_genomes/canFam3.fa",
    output:
        bed="workflow/data/canFam3.bin_200kb_all.bed",
    log:
        "workflow/data/log/canFam3_bin_bed.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        bash workflow/scripts/utils/generate_bin_bed.sh {input.fasta} {output.bed} 200000 > {log} 2>&1
        """


rule generate_canFam4_bin_bed:
    input:
        fasta="workflow/data/ref_genomes/canFam4.fa",
    output:
        bed="workflow/data/canFam4.bin_200kb_all.bed",
    log:
        "workflow/data/log/canFam4_bin_bed.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        bash workflow/scripts/utils/generate_bin_bed.sh {input.fasta} {output.bed} 200000 > {log} 2>&1
        """


rule generate_canFam3_gc_matrix:
    input:
        fasta="workflow/data/ref_genomes/canFam3.fa",
        bin_bed="workflow/data/canFam3.bin_200kb_all.bed",
    output:
        gc_matrix="workflow/data/GC/canFam3.GC_matrix.txt.gz",
    log:
        "workflow/data/log/canFam3_gc_matrix.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        python workflow/scripts/utils/generate_gc_matrix.py {input.bin_bed} {input.fasta} {output.gc_matrix} > {log} 2>&1
        """


rule generate_canFam4_gc_matrix:
    input:
        fasta="workflow/data/ref_genomes/canFam4.fa",
        bin_bed="workflow/data/canFam4.bin_200kb_all.bed",
    output:
        gc_matrix="workflow/data/GC/canFam4.GC_matrix.txt.gz",
    log:
        "workflow/data/log/canFam4_gc_matrix.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        python workflow/scripts/utils/generate_gc_matrix.py {input.bin_bed} {input.fasta} {output.gc_matrix} > {log} 2>&1
        """


rule download_T2T_tarball:
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
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
    localrule: True
    input:
        # HTTP.remote(
        storage.http(
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


# rule download_scnova_data:
#     input:
#         ancient(
#             # HTTP.remote(
#             storage.http(
#                 "https://zenodo.org/record/7697400/files/scNOVA_data_models.zip",
#                 keep_local=True,
#             )
#         ),
#     output:
#         "workflow/data/scNOVA/utils/bin_chr_length.bed",
#         "workflow/data/scNOVA/utils/bin_Genebody_all.bed",
#         "workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort_ann_sort_GC_ensemble.txt",
#         "workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort.txt",
#         "workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
#         "workflow/data/scNOVA/utils/bin_Genes_for_CNN_sort.txt.corrected",
#         "workflow/data/scNOVA/utils/Deeptool_Genes_for_CNN_merge_sort_lab_final.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_CpG_orientation_impute.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_CpG_orientation.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_GC_orientation_impute.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_GC_orientation.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_RT_orientation.txt",
#         "workflow/data/scNOVA/utils/Features_reshape_size_orientation.txt",
#         "workflow/data/scNOVA/utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
#         "workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
#         "workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort.bed",
#         "workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
#         "workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SVM.txt",
#         "workflow/data/scNOVA/utils/Strand_seq_matrix_TES_for_SVM.txt",
#         "workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
#     log:
#         touch("log/config/dl_arbigent_mappability_track.ok"),
#     conda:
#         "../envs/scNOVA/scNOVA_DL.yaml"
#     # container:
#     #     None
#     shell:
#         """
#         directory="workflow/data/ref_genomes/"
#         mkdir -p "$directory"
#         mv {input} workflow/data/scNOVA_data_models.zip
#         unzip workflow/data/scNOVA_data_models.zip -d workflow/data/
#         """
