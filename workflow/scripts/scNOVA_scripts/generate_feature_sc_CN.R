args=commandArgs(trailingOnly=TRUE)

Deeptool_result_final <- read.table(args[1], header=TRUE, sep ='\t', comment.char = "")
CNN_matrix <- Deeptool_result_final[,1:4]
DHS_matrix_for_CN <- cbind(c(1:nrow(CNN_matrix)), CNN_matrix)


#1) Calculate copy number for CN for 150 bins of 19770 genes

sv_calls <- read.table(args[2], sep = '\t', header=T)

sv_calls$X <- NA
sv_calls$CN <- matrix(2, nrow(sv_calls), 1)
sv_calls$CN[sv_calls$sv_call_name=="del_h1",1] <- 1
sv_calls$CN[sv_calls$sv_call_name=="del_h2",1] <- 1
sv_calls$CN[sv_calls$sv_call_name=="del_hom",1] <- 0
sv_calls$CN[sv_calls$sv_call_name=="dup_h1",1] <- 3
sv_calls$CN[sv_calls$sv_call_name=="dup_h2",1] <- 3
sv_calls$CN[sv_calls$sv_call_name=="dup_hom",1] <- 4
sv_calls$CN[sv_calls$sv_call_name=="idup_h1",1] <- 3
sv_calls$CN[sv_calls$sv_call_name=="idup_h2",1] <- 3

cell_list <- sort(unique(sv_calls$cell))
SC_CN_DHS <- matrix(2, nrow(DHS_matrix_for_CN), length(cell_list))


for (k in 1:length(cell_list)){


	data1 <- sv_calls[sv_calls$cell==cell_list[k],]
	data1$color <- 0
	ash12rainbow <- c("#77AADD", "#77AADD", "#114477", "#CC99BB", "#CC99BB", "#771155", "#DDDD77", "#DDDD77", "#777711", "#DDAA77", "#DDAA77", "#774411")#, "#C1FFC1")
	sv_call_name <- c("del_h1",  "del_h2",  "del_hom", "dup_h1",  "dup_h2",  "dup_hom", "inv_h1",  "inv_h2",  "inv_hom", "idup_h1", "idup_h2", "complex")#, "imputed")
	for (j in 1:length(sv_call_name)){
		tmp <- which(data1[,9]==sv_call_name[j])
		data1[tmp,15] <- j
	}

	data1_pos <- data1[,1:3]
	data1_pos_uniq <- unique(data1_pos)
	data1_pos_uniq_sort <- data1_pos_uniq[data1_pos_uniq$chrom=="chr1",]
	chrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

	for (i in 2:length(chrom)){
		data1_pos_uniq_sort<-rbind(data1_pos_uniq_sort, data1_pos_uniq[data1_pos_uniq$chrom==chrom[i],])
	}
	data1_pos_uniq_sort$posind <- c(1:nrow(data1_pos_uniq_sort))

	data1_cell <- data1[,5]
	data1_cell_uniq <- unique(data1_cell)
	data1_cell_uniq_sort <- as.matrix(sort(data1_cell_uniq))

	result <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
	result_sv <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
	result_cn <- matrix(2, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))

	for (i in 1:nrow(data1)){
		pos_ind <- which(data1_pos_uniq_sort[,1]==data1[i,1] & data1_pos_uniq_sort[,2]==data1[i,2] & data1_pos_uniq_sort[,3]==data1[i,3])
		cell_ind <- which(data1_cell_uniq_sort[,1]==data1[i,5])
		result[cell_ind, pos_ind] <- data1[i,13]
		result_sv[cell_ind, pos_ind] <- data1[i,15]
		result_cn[cell_ind, pos_ind] <- data1[i,17]
	}
	rownames(result) <- data1_cell_uniq_sort
	colnames(result) <- data1_pos_uniq_sort$posind
	result[result==Inf] <- max(result[result!=Inf])
	rownames(result_sv) <- data1_cell_uniq_sort
	colnames(result_sv) <- data1_pos_uniq_sort$posind
	rownames(result_cn) <- data1_cell_uniq_sort
	colnames(result_cn) <- data1_pos_uniq_sort$posind

	##Imputation (for single-cell CN, this will not be applied)
	imput_cand <- matrix(0, 1, 3)
	imput_cand <- imput_cand[-1,]

	#sv_call_name <- c("del_h1",  "del_h2",  "del_hom", "dup_h1",  "dup_h2",  "dup_hom", "inv_h1",  "inv_h2",  "inv_hom", "idup_h1", "idup_h2", "complex")#, "imputed")
	result_sv[result_sv==2] <- 1 #del_het (CN==1)
	result_sv[result_sv==5] <- 4 #dup_het (CN==3)
	result_sv[result_sv==10] <- 4 #idup_het (CN==3)
	result_sv[result_sv==11] <- 4 #idup_het (CN==3)
	#result_sv[result_sv==11] <- 10 #idup_het (CN==3)

	data1_pos_uniq_sort$CN <- 2
	for (i in 1:nrow(data1_pos_uniq_sort)){
		#if ((sum(result_sv[,i]==1)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 1}
		#if ((sum(result_sv[,i]==3)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 0}
		#if ((sum(result_sv[,i]==4)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 3}
		#if ((sum(result_sv[,i]==6)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 4}
		#if ((sum(result_sv[,i]==10)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 3}

		#Calculate allele frequency without complex genotype cells
		if ((sum(result_sv[,i]!=12 & result_cn[,i]==0)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 0}
		if ((sum(result_sv[,i]!=12 & result_cn[,i]==1)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 1}
		if ((sum(result_sv[,i]!=12 & result_cn[,i]==3)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 3}
		if ((sum(result_sv[,i]!=12 & result_cn[,i]==4)/nrow(result_sv))>0.8){data1_pos_uniq_sort[i,5] <- 4}
	}

	data1_pos_uniq_sort_imput <- as.data.frame(matrix(0, 1, 5))
	colnames(data1_pos_uniq_sort_imput) <- colnames(data1_pos_uniq_sort)
	for (i in 1:length(chrom)){
		if (sum(imput_cand[,1]==i)==0){data1_pos_uniq_sort_imput <- rbind(data1_pos_uniq_sort_imput, data1_pos_uniq_sort[data1_pos_uniq_sort[,1]==chrom[i],])}
		if (sum(imput_cand[,1]==i)>0){
			tmp <- data1_pos_uniq_sort[data1_pos_uniq_sort[,1]==chrom[i],]
			tmp_imput <- imput_cand[imput_cand[,1]==i,]
			for (j in 1:nrow(tmp)){
				if(j<tmp_imput[2] | j>tmp_imput[3]){data1_pos_uniq_sort_imput <- rbind(data1_pos_uniq_sort_imput, tmp[j,])}
				if(j==tmp_imput[2]){
					tmp_imput_data <- tmp[tmp_imput[2]:tmp_imput[3],]
					tmp_imput_line <- as.data.frame(cbind(as.character(tmp_imput_data[1,1]), min(tmp_imput_data[,2]), max(tmp_imput_data[,3]), min(tmp_imput_data[,4]), max(tmp_imput_data[,5])))
					colnames(tmp_imput_line) <- colnames(data1_pos_uniq_sort)
					data1_pos_uniq_sort_imput <- rbind(data1_pos_uniq_sort_imput, tmp_imput_line)
				}
			}
		}
	}
	data1_pos_uniq_sort_imput <- data1_pos_uniq_sort_imput[-1,]
	data1_pos_uniq_sort_imput[,2] <- as.integer(data1_pos_uniq_sort_imput[,2])
	data1_pos_uniq_sort_imput[,3] <- as.integer(data1_pos_uniq_sort_imput[,3])
	data1_pos_uniq_sort_imput[,4] <- as.numeric(data1_pos_uniq_sort_imput[,4])
	data1_pos_uniq_sort_imput[,5] <- as.numeric(data1_pos_uniq_sort_imput[,5])

	data1_pos_uniq_sort_imput_manual <- data1_pos_uniq_sort_imput
	data1_pos_uniq_sort_imput_manual[data1_pos_uniq_sort_imput[,5]==0,5] <- 2
	data1_pos_uniq_sort_imput <- data1_pos_uniq_sort_imput_manual


	for (j in 1:nrow(DHS_matrix_for_CN)){
		tmp <- data1_pos_uniq_sort_imput[data1_pos_uniq_sort_imput[,1]==as.character(DHS_matrix_for_CN[j,2]),]
		CN_DHS <- rbind(tmp[tmp[,2]>DHS_matrix_for_CN[j,3] & tmp[,3]<DHS_matrix_for_CN[j,4],],
		tmp[tmp[,2]<DHS_matrix_for_CN[j,3] & tmp[,3]>DHS_matrix_for_CN[j,3],],
		tmp[tmp[,2]<DHS_matrix_for_CN[j,4] & tmp[,3]>DHS_matrix_for_CN[j,4],])
		if (nrow(CN_DHS)>0){
			SC_CN_DHS[j,k] <- mean(CN_DHS[,5])
		}
		#cat(paste0(j, ' '))
	}

	library(pracma)
	CN_matrix_CNN_reshape <- t(Reshape(SC_CN_DHS[,k], 150, (19770-13)))

	CNN_features_annot <- read.table(args[3], header=T, sep ='\t', comment.char = "")
	table_original <- CN_matrix_CNN_reshape
	table_final <- as.data.frame(matrix(0, nrow(table_original), ncol(table_original)))
	for (i in 1:nrow(table_original)){
		if (CNN_features_annot[i,6]=="+"){table_final[i,] <- table_original[i,]}
		if (CNN_features_annot[i,6]=="-"){
			tmp <- table_original[i,c(150:1)]
			colnames(tmp) <- colnames(table_original)
			table_final[i,] <- tmp
			}
		#cat(paste0(i, ' '))
	}
	outname <- paste0("result_features_sc_CN/Features_reshape_", cell_list[k], "_SC_orientation_CN.txt")
	write.table(table_final, outname, row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
	# cat(paste0(k, ' '))
}
write.table(SC_CN_DHS, args[4], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
