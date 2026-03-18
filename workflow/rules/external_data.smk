# Get reference base directory from config (supports multi-user HPC setups)
REF_BASE_DIR = config.get("reference_base_dir", "workflow/data/ref_genomes")


rule dl_example_data:
    localrule: True
    output:
        touch("config/dl_example_data.ok"),
    log:
        touch("log/config/dl_example_data.ok"),
    params:
        url="https://sandbox.zenodo.org/record/1074721/files/TEST_EXAMPLE_DATA.zip",
    shell:
        """
        wget -q -O TEST_EXAMPLE_DATA.zip {params.url}
        unzip TEST_EXAMPLE_DATA.zip -d .
        """


rule download_hg19_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/hg19.fa",
    log:
        f"{REF_BASE_DIR}/log/hg19.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule download_hg38_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/hg38.fa",
    log:
        f"{REF_BASE_DIR}/log/hg38.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule download_T2T_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/T2T.fa",
    log:
        f"{REF_BASE_DIR}/log/T2T.ok",
    params:
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule download_mm10_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/mm10.fa",
    log:
        f"{REF_BASE_DIR}/log/mm10.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule download_mm39_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/mm39.fa",
    log:
        f"{REF_BASE_DIR}/log/mm39.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


# ============================================================
# Download pre-built BWA indexes from iGenomes (fast)
# ============================================================


if config.get("download_prebuilt_indexes", True):

    _igenomes_base = config["references_data"][config["reference"]].get("igenomes_base")
    assert _igenomes_base is not None, (
        f"Cannot use download_prebuilt_indexes=True with reference genome '{config['reference']}': "
        f"no iGenomes pre-built indexes are available for this genome. "
        f"Set download_prebuilt_indexes=False to build indexes locally."
    )

    rule ashleys_download_bwa_indexes:
        """Download pre-built BWA indexes directly from iGenomes to output location."""
        localrule: True
        output:
            amb=f"{get_reference_fasta()}.amb",
            ann=f"{get_reference_fasta()}.ann",
            bwt=f"{get_reference_fasta()}.bwt",
            pac=f"{get_reference_fasta()}.pac",
            sa=f"{get_reference_fasta()}.sa",
        log:
            f"{get_reference_fasta()}.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            igenomes_base=_igenomes_base,
        shell:
            """
            mkdir -p $(dirname {output.amb})
            wget -q -O {output.amb} {params.igenomes_base}/BWAIndex/genome.fa.amb
            wget -q -O {output.ann} {params.igenomes_base}/BWAIndex/genome.fa.ann
            wget -q -O {output.bwt} {params.igenomes_base}/BWAIndex/genome.fa.bwt
            wget -q -O {output.pac} {params.igenomes_base}/BWAIndex/genome.fa.pac
            wget -q -O {output.sa}  {params.igenomes_base}/BWAIndex/genome.fa.sa
            """

    rule ashleys_download_faidx:
        """Download samtools faidx directly from iGenomes to output location."""
        localrule: True
        output:
            f"{get_reference_fasta()}.fai",
        log:
            f"{get_reference_fasta()}.fai.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            igenomes_base=_igenomes_base,
        shell:
            """
            mkdir -p $(dirname {output})
            wget -q -O {output} {params.igenomes_base}/WholeGenomeFasta/genome.fa.fai
            """


rule download_canFam3_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/canFam3.fa",
    log:
        f"{REF_BASE_DIR}/log/canFam3.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule download_canFam4_reference:
    localrule: True
    output:
        f"{REF_BASE_DIR}/canFam4.fa",
    log:
        f"{REF_BASE_DIR}/log/canFam4.ok",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.fa.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output}.gz {params.url}
        gunzip {output}.gz
        touch {log}
        """


rule generate_canFam3_bin_bed:
    input:
        fasta=f"{REF_BASE_DIR}/canFam3.fa",
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
        fasta=f"{REF_BASE_DIR}/canFam4.fa",
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
        fasta=f"{REF_BASE_DIR}/canFam3.fa",
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
        fasta=f"{REF_BASE_DIR}/canFam4.fa",
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
    output:
        f"{REF_BASE_DIR}/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
    log:
        f"{REF_BASE_DIR}/log/T2T_tarball.ok",
    params:
        url="https://zenodo.org/record/7697400/files/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})
        wget -q -O {output} {params.url}
        touch {log}
        """


rule download_arbigent_mappability_track:
    localrule: True
    output:
        config["arbigent_data"]["arbigent_mapability_track"],
    log:
        touch("log/config/dl_arbigent_mappability_track.ok"),
    params:
        url="https://zenodo.org/record/7697400/files/mapping_counts_allchrs_hg38.txt",
    shell:
        """
        mkdir -p $(dirname {output})
        wget -q -O {output} {params.url}
        """


# rule download_scnova_data:
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
#     shell:
#         """
#         wget -q -O scNOVA_data_models.zip https://zenodo.org/record/7697400/files/scNOVA_data_models.zip
#         unzip scNOVA_data_models.zip -d workflow/data/
#         """
