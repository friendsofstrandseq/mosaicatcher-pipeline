args=commandArgs(trailingOnly=TRUE)

TSS_matrix <- read.table(args[1], header=F, sep ='\t', comment.char = "")
Pred_result <- read.table(args[2], header=T, sep =',', comment.char = "")
chromosome <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')

TSS_matrix_all <- data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("chr", "start", "end", "GeneID", "a", "b"))), stringsAsFactors=F)
for (i in 1:length(chromosome)){
	TSS_matrix_chr <- TSS_matrix[TSS_matrix[,1] == chromosome[i],]
	TSS_matrix_all <- rbind(TSS_matrix_all, TSS_matrix_chr)
}

Pred_all <- cbind(TSS_matrix_all, Pred_result[,2:3])
colnames(Pred_all) <- c("chr", "start", "end", "GeneID", "a", "b", "Unexpressed", "Expressed")
write.table(Pred_all, args[6], row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)


TSS_matrix <- read.table(args[1], header=F, sep ='\t', comment.char = "")
Pred_result <- read.table(args[3], header=T, sep =',', comment.char = "")
chromosome <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')

TSS_matrix_all <- data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("chr", "start", "end", "GeneID", "a", "b"))), stringsAsFactors=F)
for (i in 1:length(chromosome)){
	TSS_matrix_chr <- TSS_matrix[TSS_matrix[,1] == chromosome[i],]
	TSS_matrix_all <- rbind(TSS_matrix_all, TSS_matrix_chr)
}

Pred_all <- cbind(TSS_matrix_all, Pred_result[,2:3])
colnames(Pred_all) <- c("chr", "start", "end", "GeneID", "a", "b", "Unexpressed", "Expressed")
write.table(Pred_all, args[7], row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)



TSS_matrix <- read.table(args[1], header=F, sep ='\t', comment.char = "")
Pred_result <- read.table(args[4], header=T, sep =',', comment.char = "")
chromosome <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')

TSS_matrix_all <- data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("chr", "start", "end", "GeneID", "a", "b"))), stringsAsFactors=F)
for (i in 1:length(chromosome)){
	TSS_matrix_chr <- TSS_matrix[TSS_matrix[,1] == chromosome[i],]
	TSS_matrix_all <- rbind(TSS_matrix_all, TSS_matrix_chr)
}

Pred_all <- cbind(TSS_matrix_all, Pred_result[,2:3])
colnames(Pred_all) <- c("chr", "start", "end", "GeneID", "a", "b", "Unexpressed", "Expressed")
write.table(Pred_all, args[8], row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)


TSS_matrix <- read.table(args[1], header=F, sep ='\t', comment.char = "")
Pred_result <- read.table(args[5], header=T, sep =',', comment.char = "")
chromosome <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX')

TSS_matrix_all <- data.frame(matrix(vector(), 0, 6, dimnames=list(c(), c("chr", "start", "end", "GeneID", "a", "b"))), stringsAsFactors=F)
for (i in 1:length(chromosome)){
	TSS_matrix_chr <- TSS_matrix[TSS_matrix[,1] == chromosome[i],]
	TSS_matrix_all <- rbind(TSS_matrix_all, TSS_matrix_chr)
}

Pred_all <- cbind(TSS_matrix_all, Pred_result[,2:3])
colnames(Pred_all) <- c("chr", "start", "end", "GeneID", "a", "b", "Unexpressed", "Expressed")
write.table(Pred_all, args[9], row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)

