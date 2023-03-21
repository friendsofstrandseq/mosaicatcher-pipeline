args=commandArgs(trailingOnly=TRUE)

## 1) Load feature sets
##-------------------------------------------------------------------------------------
##Features: Sequence : GC%, CpG%, RT
##-------------------------------------------------------------------------------------

TSS_matrix <- read.table(args[1], header=TRUE, sep ='\t')
table_GC_imput <- read.table(args[2], header=F, sep ='\t', comment.char = "")
table_CpG_imput <- read.table(args[3], header=F, sep ='\t', comment.char = "")
table_RT <- read.table(args[4], header=F, sep ='\t', comment.char = "")


##-------------------------------------------------------------------------------------
##Features: Nucleosome occupancy 150 bins and copy-number normalization
##-------------------------------------------------------------------------------------
table_mononuc_norm_data1 <- read.table(args[5], header=F, sep ='\t', comment.char = "")

##Normalization by copy number (19757 X 150 copy number matrix)
CN_result_data1 <- read.table(args[6], sep = '\t', header=F)

#table_mononuc_norm_data1_cn <- table_mononuc_norm_data1/CN_result_data1 ##This was changed for ploidy mode
table_mononuc_norm_data1_cn <- table_mononuc_norm_data1/2
##-------------------------------------------------------------------------------------
##Features: Nucleosome occupancy Residual of the CV square 150 bins
##-------------------------------------------------------------------------------------

table_mononuc_var_data1 <- read.table(args[7], header=F, sep ='\t', comment.char = "")




## 2) Load gene expression information (Expressed / Unexpressed)

FPKM <- read.table(args[8], header=T, sep ='\t', comment.char = "")
TSS_matrix_woM <- TSS_matrix[TSS_matrix[,2]!="chrM",]
FPKM_woM <- FPKM[TSS_matrix[,2]!="chrM",]
FPKM_woM_LCL <- cbind(as.matrix(rowMeans(FPKM_woM[,1:9])))
Expression_label <- matrix(0, nrow(FPKM_woM_LCL), 1)
Expression_label[FPKM_woM_LCL[,1]>1,1] <- 1



## 3) Generate feature sets and target vector (Expressed = 1 / Unexpressed = 0)


Features_label <- as.matrix(Expression_label[,1])
TSS_matrix_woM_all <- TSS_matrix_woM
table_RT_all <- table_RT



##Alternative way to make input format for the transposed features (150X3)
Features_t1 <- cbind(table_mononuc_norm_data1_cn[,1], table_mononuc_var_data1[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
for (i in 2:150){
	Features_t1 <- cbind(Features_t1, cbind(table_mononuc_norm_data1_cn[,i], table_mononuc_var_data1[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
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


##This is to practice leave one chromosome validation
write.table(standard_svm_RT[,1:750], args[9], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
write.table(standard_svm_RT[,751], args[10], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
write.table(TSS_matrix_woM_all_RT[,c(2:5, 60:61)], args[11], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


