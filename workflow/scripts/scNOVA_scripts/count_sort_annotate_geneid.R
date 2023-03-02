args=commandArgs(trailingOnly=TRUE)

NO_table <- read.table(args[1], header=T, sep ='\t', comment.char = "")
GB_matrix <- read.table(args[2], header=T, sep ='\t', comment.char = "")
NO_table_annot <- cbind(NO_table[,1:3], GB_matrix$name, NO_table[,4:ncol(NO_table)])

write.table(NO_table_annot, args[3], row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)

