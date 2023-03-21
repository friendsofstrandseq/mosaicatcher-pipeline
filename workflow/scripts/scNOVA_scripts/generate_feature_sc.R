args=commandArgs(trailingOnly=TRUE)

output_filename <- args[9]

pdf(output_filename, width = 11, height = 10)

prefix <- strsplit(output_filename, "scNOVA_input_user")[[1]][1]
prefix <- substring(prefix, 1, nchar(prefix) - 1)
if (nchar(prefix) == 0) {
    prefix <- "."
}
print(prefix)

library(pracma)
Deeptool_result_final <- read.table(args[1], header=TRUE, sep ='\t', comment.char = "")
CNN_features_annot <- read.table(args[2], header=T, sep ='\t', comment.char = "")
table_CpG <- read.table(args[3], header=F, sep ='\t', comment.char = "")
table_GC <- read.table(args[4], header=F, sep ='\t', comment.char = "")
table_size <- read.table(args[5], header=F, sep ='\t', comment.char = "")  
TSS_matrix <- read.table(args[6], header=TRUE, sep ='\t')
Deeptool_mapped <- read.table(args[7], header=TRUE, sep ='\t', comment.char = "")





##Sort the order of single-cells
Deeptool_result_final_new <- Deeptool_result_final[,5:ncol(Deeptool_result_final)]
Deeptool_mapped_new <- Deeptool_mapped[,4:ncol(Deeptool_mapped)]


subclone_list <- read.table("scNOVA_input_user/input_subclonality.txt", header=T, sep ='\t', comment.char = "")
subclone_list <- read.table("scNOVA_input_user/input_subclonality.txt", header=T, sep ='\t', comment.char = "")

Deeptool_result_name <- as.data.frame(as.matrix(colnames(Deeptool_result_final_new)))
Deeptool_result_name$index <- 0
for (j in 1:nrow(Deeptool_result_name)){
	Deeptool_result_name[j,1] <- strsplit(Deeptool_result_name[j,1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
	Deeptool_result_name[j,2] <- which(subclone_list[,1]==Deeptool_result_name[j,1])
}
Deeptool_result_final_new <- Deeptool_result_final_new[,order(Deeptool_result_name[,2])]
Deeptool_result_final <- cbind(Deeptool_result_final[,1:4], Deeptool_result_final_new)


Deeptool_mapped_name <- as.data.frame(as.matrix(colnames(Deeptool_mapped_new)))
Deeptool_mapped_name$index <- 0
for (j in 1:nrow(Deeptool_mapped_name)){
	Deeptool_mapped_name[j,1] <- strsplit(Deeptool_mapped_name[j,1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
	Deeptool_mapped_name[j,2] <- which(subclone_list[,1]==Deeptool_mapped_name[j,1])
}
Deeptool_mapped_new <- Deeptool_mapped_new[,order(Deeptool_mapped_name[,2])]
Deeptool_mapped <- cbind(Deeptool_mapped[,1:3], Deeptool_mapped_new)


#This is only needed for the plotting script
FPKM <- read.table(args[8], header=T, sep ='\t', comment.char = "")
TSS_matrix_woM <- TSS_matrix[TSS_matrix[,2]!="chrM",]
FPKM_woM <- FPKM[TSS_matrix[,2]!="chrM",]

table_label <- as.matrix(colnames(Deeptool_result_final))
for (k in 1:(ncol(Deeptool_result_final)-4)){
	filename <- strsplit(table_label[k+4,1], '.sort.mdup.bam.sc_pre_mono_sort_for_mark_uniq.bam')[[1]][1]
	table_original <- t(Reshape(Deeptool_result_final[,(k+4)], 150, (19770-13)))
	table_final <- as.data.frame(matrix(0, nrow(table_original), ncol(table_original)))

	for (i in 1:nrow(table_original)){
		if (CNN_features_annot[i,6]=="+"){table_final[i,] <- table_original[i,]}
		if (CNN_features_annot[i,6]=="-"){
			tmp <- table_original[i,c(150:1)]
			colnames(tmp) <- colnames(table_original)
			table_final[i,] <- tmp
		}
	# cat(paste0(i, ' '))
	}


	#normalization
	table_mononuc <- table_final
	mapped_read <- sum(Deeptool_mapped[,(k+3)])

	##smooth spline fit for the length normalization

	table_size_flat <- matrix(0, 1, 1)
	for (i in 1:ncol(table_size)){
		table_size_flat <- rbind(table_size_flat, as.matrix(table_size[,i]))
	}
	table_size_flat <- table_size_flat[-1,]

	table_mononuc_flat <- matrix(0, 1, 1)
	for (i in 1:ncol(table_mononuc)){
		table_mononuc_flat <- rbind(table_mononuc_flat, as.matrix(table_mononuc[,i]))
	}
	table_mononuc_flat <- table_mononuc_flat[-1,]




		##Fit smooth spline (length vs. ratio to the mean) , alternative mode (log2 scale)
		Input_norm <- table_mononuc_flat
		Input_norm_sizenorm_mode2 <- matrix(0, length(Input_norm), 1)
		Input_norm_log <- log2(Input_norm+1) #log2 conversion with pseudocount 1


		plot(table_size_flat, (Input_norm_log-(mean(Input_norm_log))), xlab="size (bp)", ylab="log2(RPM/mean)")
		fit1<-smooth.spline(table_size_flat, (Input_norm_log-(mean(Input_norm_log))), df=16)
		lines(fit1,col="red",lwd=2)

		Input_norm_log_sub<-Input_norm_log
		Input_size_sub<-table_size_flat
		ratio_norm_log<-matrix(0, length(Input_norm_log_sub),1)
		a<-predict(fit1, Input_size_sub)
		ratio_norm_log[,1] <- Input_norm_log_sub-mean(Input_norm_log)-a$y
		Input_norm_sizenorm_mode2[,1] <- ratio_norm_log[,1] + mean(Input_norm_log)

		plot(table_size_flat, ratio_norm_log, main="size normalized", xlab="size (bp)", ylab="log2(RPM/mean)")
		fit1_norm<-smooth.spline(table_size_flat, ratio_norm_log,df=16)
		lines(fit1_norm,col="red",lwd=2)

		library(pracma)
		table_mononuc_norm_length <- Reshape(2^Input_norm_sizenorm_mode2-1, (19770-13), 150)
		table_mononuc_norm <- table_mononuc

		##plotting script
		##library size normalization to calculate RPM (read per million mapped read)

		table_mononuc_norm <- table_mononuc_norm_length
		plot(c(1:150), colMeans(table_mononuc_norm), type="o")
		plot(c(1:150), colMeans(table_mononuc_norm*1000000/mapped_read), type="o")
		table_mononuc_norm2 <- table_mononuc_norm*1000000/mapped_read

		data_exp <- rowMeans(FPKM_woM[,1:9]) ##This parameter need to be changed for each samples
		avg_multi <- matrix(0, 5, 150)
		avg_multi[1,] <-colMeans(table_mononuc_norm2[data_exp==0,])
		avg_multi[2,] <-colMeans(table_mononuc_norm2[data_exp>0 & data_exp<0.1,])
		avg_multi[3,] <-colMeans(table_mononuc_norm2[data_exp>0.1 & data_exp<=1,])
		avg_multi[4,] <-colMeans(table_mononuc_norm2[data_exp>1 & data_exp<=3,])
		avg_multi[5,] <-colMeans(table_mononuc_norm2[data_exp>3,])

		position <- 1:150

		plot(position, avg_multi[1,], type="o",col='gray', cex=0.1,lwd=3, ylim=range(avg_multi), xlab='genebody', ylab='Read count per million mapped read', main=filename)
		par(new=T)
		plot(position, avg_multi[2,], type="o",col='blue', cex=0.1,lwd=3, ylim=range(avg_multi), xlab='genebody', ylab='Read count per million mapped read', main=filename)
		par(new=T)
		plot(position, avg_multi[3,], type="o",col='green', cex=0.1,lwd=3, ylim=range(avg_multi), xlab='genebody', ylab='Read count per million mapped read', main=filename)
		par(new=T)
		plot(position, avg_multi[4,], type="o",col='orange', cex=0.1,lwd=3, ylim=range(avg_multi), xlab='genebody', ylab='Read count per million mapped read', main=filename)
		par(new=T)
		plot(position, avg_multi[5,], type="o",col='red', cex=0.1,lwd=3, ylim=range(avg_multi), xlab='genebody', ylab='Read count per million mapped read', main=filename)


##save the normalized result
write.table(table_mononuc_norm2, paste0("result_features_sc/Features_reshape_", filename, "_orientation_norm.txt"), row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)


}


dev.off()
