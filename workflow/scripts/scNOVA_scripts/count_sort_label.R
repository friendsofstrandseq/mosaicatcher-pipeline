args=commandArgs(trailingOnly=TRUE)

Deeptool_result <- read.table(args[1], header=TRUE, sep ='\t', comment.char = "")
Deeptool_result_new <- Deeptool_result[,4:ncol(Deeptool_result)]
Ref_bed <- read.table(args[2], header=F, sep ='\t', comment.char = "")

if (ncol(Deeptool_result)==4){Deeptool_result_new <- as.matrix(Deeptool_result_new)}

Deeptool_result_woM <- Deeptool_result[Deeptool_result[,1]!="chrM",]
Deeptool_result_new_woM <- Deeptool_result_new[Deeptool_result[,1]!="chrM",]
Ref_bed_woM <- Ref_bed[Ref_bed[,1]!="chrM",]

Deeptool_result_lab <- cbind(Deeptool_result_woM[,1:3], Ref_bed_woM[,4], Deeptool_result_new_woM)
write.table(Deeptool_result_lab, args[3], row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
