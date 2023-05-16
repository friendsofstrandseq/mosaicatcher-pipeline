rule filter_sv_calls:
    log:
        "{folder}/{sample}/log/filter_sv_calls/{sample}.log",
    input:
        "{folder}/{sample}/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv",
    output:
        "{folder}/{sample}/scNOVA_input_user/sv_calls.tsv",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/scNOVA_scripts/filter_sv_calls.py"


rule scNOVA_final_results:
    container:
        None
    input:
        get_scnova_final_output,
    output:
        "{folder}/{sample}/plots/final_results/scNOVA_{sample}.txt",
    log:
        "{folder}/log/final_blank_results/{sample}.log",
    conda:
        "../envs/mc_base.yaml"
    shell:
        "touch {output}"


rule generate_CN_for_CNN:
    container:
        None
    input:
        mosaiclassifier_final_results="{folder}/{sample}/plots/final_results/{sample}.txt",
        subclone="{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
        sv_calls_all="{folder}/{sample}/scNOVA_input_user/sv_calls.tsv",
        Deeptool_result_final="workflow/data/scNOVA/utils/Deeptool_Genes_for_CNN_merge_sort_lab_final.txt",
        CNN_features_annot="workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
    output:
        sv_calls_all_print="{folder}/{sample}/scNOVA_input_user/{clone}_sv_calls_all_print.txt",
        CN_result_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_CN_correct0.txt",
    params:
        generate_CN_for_CNN=config["scNOVA_scripts"]["generate_CN_for_CNN"],
    log:
        "{folder}/{sample}/log/{clone}/generate_CN_for_CNN.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.generate_CN_for_CNN} {input.subclone} {input.sv_calls_all} {input.Deeptool_result_final} {input.CNN_features_annot} {output.sv_calls_all_print} > {log}
        """


rule generate_CN_for_chromVAR:
    container:
        None
    input:
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        TES_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TES_for_SVM.txt",
        Genebody_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SVM.txt",
        DHS_matrix_resize="workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
        subclone="{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
        sv_calls_all="{folder}/{sample}/scNOVA_input_user/sv_calls.tsv",
    output:
        sv_calls_all_print="{folder}/{sample}/scNOVA_input_user/sv_calls_all_print_CREs.txt",
    params:
        generate_CN_for_chromVAR=config["scNOVA_scripts"]["generate_CN_for_chromVAR"],
    log:
        "{folder}/{sample}/log/generate_CN_for_chromVAR.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.generate_CN_for_chromVAR} {input.TSS_matrix} {input.TES_matrix} {input.Genebody_matrix} {input.DHS_matrix_resize} {input.subclone} {input.sv_calls_all} {output.sv_calls_all_print}  > {log}
        """


rule remove_low_quality_reads:
    container:
        None
    input:
        bam="{folder}/{sample}/selected/{cell}.sort.mdup.bam",
        check="workflow/data/scNOVA/log/dl.ok",
    output:
        bam_pre="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono.bam",
        bam_header="{folder}/{sample}/scNOVA_bam_modified/{cell}.header_test.sam",
    log:
        "{folder}/{sample}/log/remove_low_quality_reads/{cell}.log",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools view -H {input} > {output.bam_header} 
        samtools view -F 2304 {input.bam} | awk -f workflow/scripts/scNOVA_scripts/awk_1st.awk | cat {output.bam_header} - | samtools view -Sb - > {output.bam_pre}    
        """


rule sort_bam:
    log:
        "{folder}/{sample}/log/sort_bam/{cell}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono.bam",
    output:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark.bam",
    threads: 2
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools sort -@ {threads} -O BAM -o {output} {input}
        """


rule index_num1:
    log:
        "{folder}/{sample}/log/index_num1/{cell}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark.bam",
    output:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark.bam.bai",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools index {input}
        """


rule remove_dup:
    log:
        "{folder}/{sample}/log/remove_dup/{cell}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark.bam",
    output:
        bam_uniq="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bam_metrix="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono.metrix_dup.txt",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        bammarkduplicates markthreads=2 I={input.bam} O={output.bam_uniq} M={output.bam_metrix} index=1 rmdup=1
        """


rule index_num2:
    log:
        "{folder}/{sample}/log/index_num2/{cell}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
    output:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools index {input}
        """


rule count_reads_split:
    log:
        "{folder}/{sample}/log/count_reads_split/{cell}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bai="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_split/{cell}.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam}  -bed workflow/data/scNOVA/utils/bin_Genebody_all.bed > {output.tab}
        """


rule count_reads_split_aggr:
    log:
        "{folder}/{sample}/log/count_reads_split_aggr.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_split/{cell}.tab",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/{sample}.tab",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_sort_by_coordinate:
    log:
        "{folder}/{sample}/log/count_sort_by_coordinate/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/{sample}.tab",
    output:
        "{folder}/{sample}/scNOVA_result/{sample}_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_sort_annotate_geneid:
    log:
        "{folder}/{sample}/log/count_sort_annotate_geneid/{sample}.log",
    container:
        None
    input:
        count_table="{folder}/{sample}/scNOVA_result/{sample}_sort.txt",
        GB_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
    output:
        "{folder}/{sample}/scNOVA_result/{sample}_sort_geneid.txt",
    params:
        count_sort_annotate_geneid=config["scNOVA_scripts"][
            "count_sort_annotate_geneid"
        ],
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        Rscript {params.count_sort_annotate_geneid} {input.count_table} {input.GB_matrix} {output}  
        """


rule filter_input_subclonality:
    log:
        "{folder}/{sample}/log/filter_input_subclonality/{clone}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
    output:
        "{folder}/{sample}/scNOVA_input_user/input_subclonality_{clone}.txt",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/scNOVA_scripts/filter_input_subclonality.py"


rule merge_bam_clones:
    log:
        "{folder}/{sample}/log/merge_bam_clones/{clone}.log",
    container:
        None
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
        bai=lambda wc: expand(
            "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
        input_subclonality="{folder}/{sample}/scNOVA_input_user/input_subclonality_{clone}.txt",
    output:
        bam="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam",
        bai="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam.bai",
        subclonality_colnames="{folder}/{sample}/scNOVA_input_user/{clone}_input_subclonality_colnames.txt",
        line="{folder}/{sample}/scNOVA_input_user/{clone}_line.txt",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        perl workflow/scripts/scNOVA_scripts/merge_bam_clones.pl {input.input_subclonality} {output.subclonality_colnames} {output.line}
        """


rule count_reads_for_DNN:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN/{clone}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam",
        bai="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_for_DNN/Deeptool_Genes_for_CNN_{clone}.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam}  -bed workflow/data/scNOVA/utils/bin_Genes_for_CNN_sort.txt.corrected > {output.tab}
        """


rule count_reads_for_DNN_aggr:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_aggr/{sample}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_for_DNN/Deeptool_Genes_for_CNN_{clone}.tab",
            clone=clones[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}.tab",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_reads_for_DNN_sc:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sc/{cell}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bai="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_for_DNN_sc/Deeptool_Genes_for_CNN_{cell}.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam}  -bed workflow/data/scNOVA/utils/bin_Genes_for_CNN_sort.txt.corrected > {output.tab}
        """


rule count_reads_for_DNN_sc_aggr:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sc_aggr/{sample}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_for_DNN_sc/Deeptool_Genes_for_CNN_{cell}.tab",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc.tab",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_reads_chr_length:
    log:
        "{folder}/{sample}/log/count_reads_chr_length/{clone}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam",
        bai="{folder}/{sample}/scNOVA_bam_merge/{clone}.merge.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_chr_length/Deeptool_chr_length_{clone}.tab",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam} -bed workflow/data/scNOVA/utils/bin_chr_length.bed > {output.tab}
        """


rule count_reads_chr_length_aggr:
    log:
        "{folder}/{sample}/log/count_reads_chr_length_aggr/{sample}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_chr_length/Deeptool_chr_length_{clone}.tab",
            clone=clones[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/Deeptool_chr_length_{sample}.tab",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_reads_chr_length_sc:
    log:
        "{folder}/{sample}/log/count_reads_chr_length_sc/{cell}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bai="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_chr_length_sc/Deeptool_chr_length_{cell}.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam}  -bed workflow/data/scNOVA/utils/bin_chr_length.bed > {output.tab}
        """


rule count_reads_chr_length_sc_aggr:
    log:
        "{folder}/{sample}/log/count_reads_chr_length_sc_aggr/{sample}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_chr_length_sc/Deeptool_chr_length_{cell}.tab",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/Deeptool_chr_length_{sample}_sc.tab",
    conda:
        "../envs/mc_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_reads_for_DNN_sort:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sort/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}.tab",
    output:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_reads_for_DNN_sort_lab:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sort_lab/{sample}.log",
    container:
        None
    input:
        count_reads_sort="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort.txt",
        Ref_bed="workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort.txt",
    output:
        count_reads_sort_label="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort_lab.txt",
    params:
        count_sort_label=config["scNOVA_scripts"]["count_sort_label"],
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.count_sort_label} {input.count_reads_sort} {input.Ref_bed} {output.count_reads_sort_label}
        """


rule count_reads_for_DNN_sort_label_sort:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sort_label_sort/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort_lab.txt",
    output:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort_lab_final.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k4,4n -t$'\t' {input} > {output}
        """


rule count_reads_for_DNN_normalization:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_normalization/{clone}.log",
    container:
        None
    input:
        count_reads_chr_length="{folder}/{sample}/scNOVA_result/Deeptool_chr_length_{sample}.tab",
        count_reads_sort_label="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sort_lab.txt",
        CNN_features_annot="workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
        table_CpG="workflow/data/scNOVA/utils/Features_reshape_CpG_orientation.txt",
        table_GC="workflow/data/scNOVA/utils/Features_reshape_GC_orientation.txt",
        table_size="workflow/data/scNOVA/utils/Features_reshape_size_orientation.txt",
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        FPKM="workflow/data/scNOVA/utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
        CN_result_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_CN_correct0.txt",
    output:
        plot="{folder}/{sample}/scNOVA_result/Features_reshape_{sample}_{clone}_orientation_norm_qc.pdf",
        table_mononuc_norm_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_norm.txt",
    params:
        count_norm=config["scNOVA_scripts"]["count_norm"],
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.count_norm} {input.count_reads_chr_length} {input.count_reads_sort_label} {input.CNN_features_annot} {input.table_CpG} {input.table_GC} {input.table_size} {input.TSS_matrix} {input.FPKM} {input.CN_result_data1} {output.plot} {output.table_mononuc_norm_data1}
        """


rule count_reads_for_DNN_sc_sort:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sc_sort/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc.tab",
    output:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_reads_for_DNN_sc_sort_lab:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sc_sort_lab/{sample}.log",
    container:
        None
    input:
        count_reads_sort="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort.txt",
        Ref_bed="workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort.txt",
    output:
        count_reads_sort_label="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort_lab.txt",
    params:
        count_sort_label=config["scNOVA_scripts"]["count_sort_label"],
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.count_sort_label} {input.count_reads_sort} {input.Ref_bed} {output.count_reads_sort_label}
        """


rule count_reads_for_DNN_sc_sort_label_sort:
    log:
        "{folder}/{sample}/log/count_reads_for_DNN_sc_sort_label_sort/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort_lab.txt",
    output:
        "{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort_lab_final.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k4,4n -t$'\t' {input} > {output}
        """


rule generate_feature_sc_var:
    log:
        "{folder}/{sample}/log/generate_feature_sc_var/{clone}.log",
    container:
        None
    input:
        subclone_list="{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
        count_reads_sc_sort="{folder}/{sample}/scNOVA_result/Deeptool_Genes_for_CNN_{sample}_sc_sort_lab_final.txt",
        Ref_bed_annot="workflow/data/scNOVA/utils/bin_Genes_for_CNN_num_sort_ann_sort_GC_ensemble.txt",
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        CNN_features_annot="workflow/data/scNOVA/utils/bin_Genes_for_CNN_reshape_annot.txt",
        FPKM="workflow/data/scNOVA/utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
        CN_result_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_CN_correct0.txt",
    output:
        plot="{folder}/{sample}/scNOVA_result/Features_reshape_{sample}_{clone}_Resid_orientation_qc.pdf",
        table_mononuc_var_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_Resid_orientation.txt",
    params:
        feature_sc_var=config["scNOVA_scripts"]["feature_sc_var"],
    log:
        "{folder}/{sample}/log/generate_feature_sc_var_{clone}.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.feature_sc_var} {input.subclone_list} {input.count_reads_sc_sort} {input.Ref_bed_annot} {input.TSS_matrix} {input.CNN_features_annot} {input.FPKM} {input.CN_result_data1} {output.plot} {output.table_mononuc_var_data1} > {log} 2>&1
       """


rule combine_features:
    log:
        "{folder}/{sample}/log/combine_features/{clone}.log",
    container:
        None
    input:
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        table_GC_imput="workflow/data/scNOVA/utils/Features_reshape_GC_orientation_impute.txt",
        table_CpG_imput="workflow/data/scNOVA/utils/Features_reshape_CpG_orientation_impute.txt",
        table_RT="workflow/data/scNOVA/utils/Features_reshape_RT_orientation.txt",
        table_mononuc_norm_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_norm.txt",
        CN_result_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_orientation_CN_correct0.txt",
        table_mononuc_var_data1="{folder}/{sample}/scNOVA_result/Features_reshape_{clone}_Resid_orientation.txt",
        FPKM="workflow/data/scNOVA/utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
    output:
        features="{folder}/{sample}/scNOVA_result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{clone}.txt",
        exp="{folder}/{sample}/scNOVA_result/Expression_all_{clone}.txt",
        TSS_annot="{folder}/{sample}/scNOVA_result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
    params:
        combine_features=config["scNOVA_scripts"]["combine_features"],
    log:
        "{folder}/{sample}/log/combine_features_{clone}.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="48:00:00",
    shell:
        """
        Rscript {params.combine_features} {input.TSS_matrix} {input.table_GC_imput} {input.table_CpG_imput} {input.table_RT} {input.table_mononuc_norm_data1} {input.CN_result_data1} {input.table_mononuc_var_data1} {input.FPKM} {output.features} {output.exp} {output.TSS_annot}
        """


rule infer_expressed_genes_split:
    log:
        "{folder}/{sample}/log/infer_expressed_genes_split/{clone}_{chrom}_{i}.log",
    container:
        None
    input:
        features="{folder}/{sample}/scNOVA_result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{clone}.txt",
        TSS_annot="{folder}/{sample}/scNOVA_result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
    output:
        train="{folder}/{sample}/scNOVA_result_CNN/{chrom}/DNN_train{i}_output_ypred_{clone}.csv",
    conda:
        "../envs/scNOVA/scNOVA_DL.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/scNOVA_scripts/Deeplearning_Nucleosome_predict_train_RPE.py"


rule gather_infer_expressed_genes_split:
    log:
        "{folder}/{sample}/log/gather_infer_expressed_genes_split/{clone}_{i}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result_CNN/{chrom}/DNN_train{i}_output_ypred_{clone}.csv",
            i=wc.i,
            clone=wc.clone,
            chrom=config["chromosomes"],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        "{folder}/{sample}/scNOVA_result_CNN/DNN_train{i}_output_ypred_{clone}.csv",
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/scNOVA_scripts/gather_infer_expr_genes_split.py"


rule aggr_models_touch:
    log:
        "{folder}/{sample}/log/aggr_models_touch/{clone}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result_CNN/DNN_train{i}_output_ypred_{clone}.csv",
            i=["5", "20", "40", "80"],
            clone=wc.clone,
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        touch("{folder}/{sample}/scNOVA_result_CNN/DNN_train_models_{clone}.ok"),


rule annot_expressed_genes:
    log:
        "{folder}/{sample}/log/annot_expressed_genes/{clone}.log",
    container:
        None
    input:
        TSS_annot="{folder}/{sample}/scNOVA_result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
        train80="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_{clone}.csv",
        train40="{folder}/{sample}/scNOVA_result_CNN/DNN_train40_output_ypred_{clone}.csv",
        train20="{folder}/{sample}/scNOVA_result_CNN/DNN_train20_output_ypred_{clone}.csv",
        train5="{folder}/{sample}/scNOVA_result_CNN/DNN_train5_output_ypred_{clone}.csv",
        check="{folder}/{sample}/scNOVA_result_CNN/DNN_train_models_{clone}.ok",
    output:
        train80_annot="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_{clone}_annot.txt",
        train40_annot="{folder}/{sample}/scNOVA_result_CNN/DNN_train40_output_ypred_{clone}_annot.txt",
        train20_annot="{folder}/{sample}/scNOVA_result_CNN/DNN_train20_output_ypred_{clone}_annot.txt",
        train5_annot="{folder}/{sample}/scNOVA_result_CNN/DNN_train5_output_ypred_{clone}_annot.txt",
    params:
        annot_expressed=config["scNOVA_scripts"]["annot_expressed"],
    log:
        "{folder}/{sample}/log/annot_expressed_genes_{clone}.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        Rscript {params.annot_expressed} {input.TSS_annot} {input.train80} {input.train40} {input.train20} {input.train5} {output.train80_annot} {output.train40_annot} {output.train20_annot} {output.train5_annot}
        """


rule infer_differential_gene_expression:
    log:
        "{folder}/{sample}/log/infer_differential_gene_expression/{sample}.log",
    container:
        None
    input:
        Genebody_NO="{folder}/{sample}/scNOVA_result/{sample}_sort.txt",
        clonality="{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        GB_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
        CNN_result1="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_clone1_annot.txt",
        CNN_result2="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_clone2_annot.txt",
        input_matrix="workflow/data/scNOVA_input_SV_affected_genes.txt",
    output:
        pdf="{folder}/{sample}/scNOVA_result_plots/Result_scNOVA_plots_{sample}.pdf",
        final_result="{folder}/{sample}/scNOVA_result_CNN/Expressed_train80_final_result.txt",
    params:
        infer_diff_gene_expression=config["scNOVA_scripts"][
            "infer_diff_gene_expression"
        ],
    log:
        "{folder}/{sample}/log/infer_diff_gene_expression.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    shell:
        """
        Rscript {params.infer_diff_gene_expression} {input.Genebody_NO} {input.clonality} {input.TSS_matrix} {input.GB_matrix} {input.CNN_result1} {input.CNN_result2} {input.input_matrix} {output.pdf} {output.final_result}
        """


rule infer_differential_gene_expression_alt:
    log:
        "{folder}/{sample}/log/infer_differential_gene_expression_alt/{sample}.log",
    container:
        None
    input:
        Genebody_NO="{folder}/{sample}/scNOVA_result/{sample}_sort.txt",
        clonality="{folder}/{sample}/scNOVA_input_user/input_subclonality.txt",
        TSS_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_TSS_for_SVM.txt",
        GB_matrix="workflow/data/scNOVA/utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
        CNN_result1="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_clone1_annot.txt",
        CNN_result2="{folder}/{sample}/scNOVA_result_CNN/DNN_train80_output_ypred_clone2_annot.txt",
        input_matrix="workflow/data/scNOVA_input_SV_affected_genes.txt",
        final_result="{folder}/{sample}/scNOVA_result_CNN/Expressed_train80_final_result.txt",
    output:
        result_table="{folder}/{sample}/scNOVA_result/result_PLSDA_{sample}.txt",
        result_plot="{folder}/{sample}/scNOVA_result_plots/Result_scNOVA_plots_{sample}_alternative_PLSDA.pdf",
    params:
        infer_diff_gene_expression_alt=config["scNOVA_scripts"][
            "infer_diff_gene_expression_alt"
        ],
    log:
        "{folder}/{sample}/log/infer_diff_gene_expression_alt.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    shell:
        """
        Rscript {params.infer_diff_gene_expression_alt} {input.Genebody_NO} {input.clonality} {input.TSS_matrix} {input.GB_matrix} {input.CNN_result1} {input.CNN_result2} {input.input_matrix} {output.result_table} {output.result_plot} {input.final_result} 
        """


rule count_reads_CREs:
    log:
        "{folder}/{sample}/log/count_reads_CREs/{cell}.log",
    container:
        None
    input:
        bam="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bai="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai",
    output:
        tab="{folder}/{sample}/scNOVA_result/count_reads_CREs/{cell}_CREs_2kb.tab",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam}  -bed workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort.bed > {output.tab}
        """


rule count_reads_CREs_aggr:
    log:
        "{folder}/{sample}/log/count_reads_CREs_aggr/{sample}.log",
    container:
        None
    input:
        lambda wc: expand(
            "{folder}/{sample}/scNOVA_result/count_reads_CREs/{cell}_CREs_2kb.tab",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        tab="{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
    conda:
        "../envs/mc_base.yaml"
    script:
        "../scripts/scNOVA_scripts/dev_aggr.py"


rule count_sort_by_coordinate_CREs:
    log:
        "{folder}/{sample}/log/count_sort_by_coordinate_CREs/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb.tab",
    output:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_sort_annotate_chrid_CREs:
    log:
        "{folder}/{sample}/log/count_sort_annotate_chrid_CREs/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort.txt",
    output:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort_num.txt",
    params:
        count_sort_annotate_chrid_CREs=config["scNOVA_scripts"][
            "count_sort_annotate_chrid_CREs"
        ],
    log:
        "{folder}/{sample}/log/count_sort_annotate_chrid_CREs.log",
    conda:
        "../envs/scNOVA/scNOVA_R.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        Rscript {params.count_sort_annotate_chrid_CREs} {input} {output} 
        """


rule count_sort_annotate_chrid_CREs_sort:
    log:
        "{folder}/{sample}/log/count_sort_annotate_chrid_CREs_sort/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort_num.txt",
    output:
        "{folder}/{sample}/scNOVA_result/{sample}_CREs_2kb_sort_num_sort_for_chromVAR.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1n -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule split_bam_WC:
    log:
        "{folder}/{sample}/log/split_bam_WC/{cell}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
    output:
        bam_header="{folder}/{sample}/scNOVA_bam_modified/{cell}.header_WC.sam",
        bam_C1="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C1.bam",
        bam_C2="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C2.bam",
        bam_W1="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W1.bam",
        bam_W2="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W2.bam",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools view -H {input} > {output.bam_header}
        samtools view -f 99 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_C1}
        samtools view -f 147 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_C2}
        samtools view -f 83 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_W1}
        samtools view -f 163 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_W2}
        """


rule split_bam_WC_merge:
    log:
        "{folder}/{sample}/log/split_bam_WC_merge/{cell}.log",
    container:
        None
    input:
        bam_C1="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C1.bam",
        bam_C2="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C2.bam",
        bam_W1="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W1.bam",
        bam_W2="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W2.bam",
    output:
        bam_C="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam",
        bam_W="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools merge {output.bam_C} {input.bam_C1} {input.bam_C2}
        samtools merge {output.bam_W} {input.bam_W1} {input.bam_W2}
        """


rule split_bam_WC_index:
    log:
        "{folder}/{sample}/log/split_bam_WC_index/{cell}.log",
    container:
        None
    input:
        bam_C="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam",
        bam_W="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam",
    output:
        bam_C_ind="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam.bai",
        bam_W_ind="{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam.bai",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        samtools index {input.bam_C}
        samtools index {input.bam_W}
        """


rule perl_split_sc:
    log:
        "{folder}/{sample}/log/perl_split_sc/{sample}.log",
    container:
        None
    input:
        strandphaser_output="{folder}/{sample}/strandphaser/strandphaser_phased_haps_merged.txt",
        bam_C_ind=lambda wc: expand(
            "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam.bai",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
        bam_W_ind=lambda wc: expand(
            "{folder}/{sample}/scNOVA_bam_modified/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam.bai",
            cell=bam_per_sample_selected[wc.sample],
            sample=wc.sample,
            folder=config["data_location"],
        ),
    output:
        nucleosome_sampleA="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleA/result.H1.bam",
        nucleosome_sampleB="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleB/result.H2.bam",
        strandphaser_output_copy="{folder}/{sample}/scNOVA_input_user/strandphaser_output_copy.txt",
    log:
        "{folder}/{sample}/log/perl_split_sc.log",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        perl workflow/scripts/scNOVA_scripts/perl_test_all_snake.pl {input.strandphaser_output} {output.nucleosome_sampleA} {output.nucleosome_sampleB} {output.strandphaser_output_copy}
        """


rule count_reads_CREs_haplo:
    log:
        "{folder}/{sample}/log/count_reads_CREs_haplo/{sample}.log",
    container:
        None
    input:
        bam1="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleA/result.H1.bam",
        bam2="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleB/result.H2.bam",
    output:
        tab="{folder}/{sample}/scNOVA_result_haplo/Deeptool_DHS_2kb_H1H2.tab",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam1} {input.bam2} -bed workflow/data/scNOVA/utils/regions_all_hg38_v2_resize_2kb_sort.bed > {output.tab}
        """


rule count_reads_CREs_haplo_sort_by_coordinate:
    log:
        "{folder}/{sample}/log/count_reads_CREs_haplo_sort_by_coordinate/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result_haplo/Deeptool_DHS_2kb_H1H2.tab",
    output:
        "{folder}/{sample}/scNOVA_result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_reads_genebody_haplo:
    log:
        "{folder}/{sample}/log/count_reads_genebody_haplo/{sample}.log",
    container:
        None
    input:
        bam1="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleA/result.H1.bam",
        bam2="{folder}/{sample}/scNOVA_nucleosomes_bam/nucleosome_sampleB/result.H2.bam",
    output:
        tab="{folder}/{sample}/scNOVA_result_haplo/Deeptool_Genebody_H1H2.tab",
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    threads: 1
    conda:
        "../envs/scNOVA/scNOVA_bioinfo_tools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam1} {input.bam2}  -bed workflow/data/scNOVA/utils/bin_Genebody_all.bed > {output.tab}
        """


rule count_reads_genebody_haplo_sort_by_coordinate_genebody:
    log:
        "{folder}/{sample}/log/count_reads_genebody_haplo_sort_by_coordinate_genebody/{sample}.log",
    container:
        None
    input:
        "{folder}/{sample}/scNOVA_result_haplo/Deeptool_Genebody_H1H2.tab",
    output:
        "{folder}/{sample}/scNOVA_result_haplo/Deeptool_Genebody_H1H2_sort.txt",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """
