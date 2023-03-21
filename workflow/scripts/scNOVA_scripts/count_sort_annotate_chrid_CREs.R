args=commandArgs(trailingOnly=TRUE)

data1_resize<-read.table(args[1], sep = '\t', header=T, comment.char = "")
data1_resize$chr <- 0
data1_resize[data1_resize[,1]=="chr1",ncol(data1_resize)] <- 1
data1_resize[data1_resize[,1]=="chr2",ncol(data1_resize)] <- 2
data1_resize[data1_resize[,1]=="chr3",ncol(data1_resize)] <- 3
data1_resize[data1_resize[,1]=="chr4",ncol(data1_resize)] <- 4
data1_resize[data1_resize[,1]=="chr5",ncol(data1_resize)] <- 5
data1_resize[data1_resize[,1]=="chr6",ncol(data1_resize)] <- 6
data1_resize[data1_resize[,1]=="chr7",ncol(data1_resize)] <- 7
data1_resize[data1_resize[,1]=="chr8",ncol(data1_resize)] <- 8
data1_resize[data1_resize[,1]=="chr9",ncol(data1_resize)] <- 9
data1_resize[data1_resize[,1]=="chr10",ncol(data1_resize)] <- 10
data1_resize[data1_resize[,1]=="chr11",ncol(data1_resize)] <- 11
data1_resize[data1_resize[,1]=="chr12",ncol(data1_resize)] <- 12
data1_resize[data1_resize[,1]=="chr13",ncol(data1_resize)] <- 13
data1_resize[data1_resize[,1]=="chr14",ncol(data1_resize)] <- 14
data1_resize[data1_resize[,1]=="chr15",ncol(data1_resize)] <- 15
data1_resize[data1_resize[,1]=="chr16",ncol(data1_resize)] <- 16
data1_resize[data1_resize[,1]=="chr17",ncol(data1_resize)] <- 17
data1_resize[data1_resize[,1]=="chr18",ncol(data1_resize)] <- 18
data1_resize[data1_resize[,1]=="chr19",ncol(data1_resize)] <- 19
data1_resize[data1_resize[,1]=="chr20",ncol(data1_resize)] <- 20
data1_resize[data1_resize[,1]=="chr21",ncol(data1_resize)] <- 21
data1_resize[data1_resize[,1]=="chr22",ncol(data1_resize)] <- 22
data1_resize[data1_resize[,1]=="chrX",ncol(data1_resize)] <- 23
data1_resize <- data1_resize[,c(ncol(data1_resize), 2:(ncol(data1_resize)-1))]
write.table(data1_resize, args[2], row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)

