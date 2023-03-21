args=commandArgs(trailingOnly=TRUE)


TSS_matrix <- read.table(args[1], header=TRUE, sep ='\t')
subclonality <- read.table(args[2], header=TRUE, sep ='\t')


##-------------------------------------------------------------------------------------
##Features: Sequence : GC%, CpG%, RT
##-------------------------------------------------------------------------------------
table_GC_imput <- read.table(args[3], header=F, sep ='\t', comment.char = "")
table_CpG_imput <- read.table(args[4], header=F, sep ='\t', comment.char = "")
table_RT <- read.table(args[5], header=F, sep ='\t', comment.char = "")

##-------------------------------------------------------------------------------------
##Features: Nucleosome occupancy 150 bins and copy-number normalization
##-------------------------------------------------------------------------------------

files <- list.files("result_features_sc/", pattern="_norm.txt$", full.names=TRUE)

for (k in 1:length(files)){
	filename <- strsplit(files[k], 'result_features_sc//Features_reshape_')[[1]][2]
	filename <- strsplit(filename, '_orientation_norm.txt')[[1]][1]
	filename <- strsplit(filename, '.sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam')[[1]][1]
	filename_CN <- paste0("result_features_sc_CN/Features_reshape_", filename, "_SC_orientation_CN.txt")

	table_mononuc_norm_data1 <- read.table(files[k], header=F, sep ='\t', comment.char = "")

	##Normalization by copy number (19757 X 150 copy number matrix)
	CN_result_data1 <- read.table(filename_CN, sep = '\t', header=F)
	table_mononuc_norm_data1_cn <- table_mononuc_norm_data1/CN_result_data1


	##-------------------------------------------------------------------------------------
	##Features: Nucleosome occupancy Residual of the CV square 150 bins
	##-------------------------------------------------------------------------------------


	FPKM <- read.table(args[6], header=T, sep ='\t', comment.char = "")
	TSS_matrix_woM <- TSS_matrix[TSS_matrix[,2]!="chrM",]
	FPKM_woM <- FPKM[TSS_matrix[,2]!="chrM",]
	FPKM_woM_LCL <- cbind(as.matrix(rowMeans(FPKM_woM[,1:9])))
	Expression_label <- matrix(0, nrow(FPKM_woM_LCL), 1)
	Expression_label[FPKM_woM_LCL[,1]>1,1] <- 1



	Features_label <- as.matrix(Expression_label[,1])#For BCLL01, use mean of LCL
	TSS_matrix_woM_all <- TSS_matrix_woM
	table_RT_all <- table_RT



	##Alternative way to make input format for the transposed features (150X4)
	Features_t1 <- cbind(table_mononuc_norm_data1_cn[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
	for (i in 2:150){
		Features_t1 <- cbind(Features_t1, cbind(table_mononuc_norm_data1_cn[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
	}
	Features_tall <- rbind(Features_t1)
	Features_both_sub <- Features_tall


	standard_svm <- data.frame(Features_both_sub)
	standard_svm$Type<-rep(0, nrow(standard_svm))
	standard_svm[Features_label=="0",ncol(Features_both_sub)+1] <- 0
	standard_svm[Features_label=="1",ncol(Features_both_sub)+1] <- 1
	standard_svm$Type <- as.factor(standard_svm$Type)

	standard_svm_RT <- standard_svm[is.na(rowSums(table_RT_all))==0,]
	TSS_matrix_woM_all_RT <- TSS_matrix_woM_all[is.na(rowSums(table_RT_all))==0,]

	#standard_svm_RT_union <- rbind(standard_svm_RT_union, standard_svm_RT)
	#TSS_matrix_woM_all_RT_union <- rbind(TSS_matrix_woM_all_RT_union, TSS_matrix_woM_all_RT)

	write.table(standard_svm_RT[,1:601], paste0("result_features_sc/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_", filename, "_wovar_exp.txt"), row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
	# cat(paste0(k, ' '))
}
write.table(TSS_matrix_woM_all_RT[,c(2:5, 60:61)], args[7], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)



