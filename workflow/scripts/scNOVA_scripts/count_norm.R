args <- commandArgs(trailingOnly = TRUE)

output_filename <- args[10]

pdf(output_filename, width = 11, height = 10)

prefix <- strsplit(output_filename, "scNOVA_result")[[1]][1]
prefix <- substring(prefix, 1, nchar(prefix) - 1)
if (nchar(prefix) == 0) {
    prefix <- "."
}
print(prefix)

library(pracma)
Deeptool_chr_length <- read.table(args[1], header = TRUE, sep = "\t", comment.char = "")
Deeptool_result_final <- read.table(args[2], header = TRUE, sep = "\t", comment.char = "")

clone_num <- ncol(Deeptool_result_final) - 4

# for (k in 1:clone_num){
clone_name <- strsplit(args[9], paste0(prefix, "/scNOVA_result/Features_reshape_"))[[1]][2]
clone_name <- strsplit(clone_name, paste0(prefix, "/scNOVA_result/_orientation_CN_correct0.txt"))[[1]][1]
k <- as.numeric(strsplit(clone_name, "clone")[[1]][2])
filename <- paste0(prefix, "/scNOVA_result/Features_reshape_", clone_name, "_orientation_norm.txt")
Deeptool_result_final_reshape <- t(Reshape(Deeptool_result_final[, (k + 4)], 150, (19770 - 13)))

#----------------------------------------------------------------------------------------------------------------
# Change orientation of the features according to the strand information of the genes
#----------------------------------------------------------------------------------------------------------------
CNN_features_annot <- read.table(args[3], header = T, sep = "\t", comment.char = "")
table_original <- Deeptool_result_final_reshape
table_final <- as.data.frame(matrix(0, nrow(table_original), ncol(table_original)))
for (i in 1:nrow(table_original)) {
    if (CNN_features_annot[i, 6] == "+") {
        table_final[i, ] <- table_original[i, ]
    }
    if (CNN_features_annot[i, 6] == "-") {
        tmp <- table_original[i, c(150:1)]
        colnames(tmp) <- colnames(table_original)
        table_final[i, ] <- tmp
    }
    # cat(paste0(i, ' '))
}

#----------------------------------------------------------------------------------------------------------------
# Length normalization and make Genebody plot based on RNA levels
#----------------------------------------------------------------------------------------------------------------
# library size : This need to be calculated from entire genome
# length : This will be calculated from bin by bin
# GC : do we need this normalization or not if we have GC as a feature

mapped_read <- sum(Deeptool_chr_length[, (k + 3)]) ## This parameter need to be changed for each samples

table_mononuc <- table_final
table_CpG <- read.table(args[4], header = F, sep = "\t", comment.char = "")
table_GC <- read.table(args[5], header = F, sep = "\t", comment.char = "")
table_size <- read.table(args[6], header = F, sep = "\t", comment.char = "")
TSS_matrix <- read.table(args[7], header = TRUE, sep = "\t")

FPKM <- read.table(args[8], header = T, sep = "\t", comment.char = "")
TSS_matrix_woM <- TSS_matrix[TSS_matrix[, 2] != "chrM", ]
FPKM_woM <- FPKM[TSS_matrix[, 2] != "chrM", ]

## smooth spline fit for the length normalization

table_size_flat <- matrix(0, 1, 1)
for (i in 1:ncol(table_size)) {
    table_size_flat <- rbind(table_size_flat, as.matrix(table_size[, i]))
}
table_size_flat <- table_size_flat[-1, ]

table_mononuc_flat <- matrix(0, 1, 1)
for (i in 1:ncol(table_mononuc)) {
    table_mononuc_flat <- rbind(table_mononuc_flat, as.matrix(table_mononuc[, i]))
}
table_mononuc_flat <- table_mononuc_flat[-1, ]

## Fit smooth spline (length vs. ratio to the mean) , alternative mode (log2 scale)
Input_norm <- table_mononuc_flat
Input_norm_sizenorm_mode2 <- matrix(0, length(Input_norm), 1)
Input_norm_log <- log2(Input_norm + 1) # log2 conversion with pseudocount 1


plot(table_size_flat, (Input_norm_log - (mean(Input_norm_log))), xlab = "size (bp)", ylab = "log2(RPM/mean)")
fit1 <- smooth.spline(table_size_flat, (Input_norm_log - (mean(Input_norm_log))), df = 16)
lines(fit1, col = "red", lwd = 2)

Input_norm_log_sub <- Input_norm_log
Input_size_sub <- table_size_flat
ratio_norm_log <- matrix(0, length(Input_norm_log_sub), 1)
a <- predict(fit1, Input_size_sub)
ratio_norm_log[, 1] <- Input_norm_log_sub - mean(Input_norm_log) - a$y
Input_norm_sizenorm_mode2[, 1] <- ratio_norm_log[, 1] + mean(Input_norm_log)

plot(table_size_flat, ratio_norm_log, main = "size normalized", xlab = "size (bp)", ylab = "log2(RPM/mean)")
fit1_norm <- smooth.spline(table_size_flat, ratio_norm_log, df = 16)
lines(fit1_norm, col = "red", lwd = 2)

library(pracma)
table_mononuc_norm_length <- Reshape(2^Input_norm_sizenorm_mode2 - 1, (19770 - 13), 150)
table_mononuc_norm <- table_mononuc

## plotting script
## library size normalization to calculate RPM (read per million mapped read)

table_mononuc_norm <- table_mononuc_norm_length
plot(c(1:150), colMeans(table_mononuc_norm), type = "o")
plot(c(1:150), colMeans(table_mononuc_norm * 1000000 / mapped_read), type = "o")
table_mononuc_norm2 <- table_mononuc_norm * 1000000 / mapped_read

data_exp <- rowMeans(FPKM_woM[, 1:9]) ## This parameter need to be changed for each samples
avg_multi <- matrix(0, 5, 150)
avg_multi[1, ] <- colMeans(table_mononuc_norm2[data_exp == 0, ])
avg_multi[2, ] <- colMeans(table_mononuc_norm2[data_exp > 0 & data_exp < 0.1, ])
avg_multi[3, ] <- colMeans(table_mononuc_norm2[data_exp > 0.1 & data_exp <= 1, ])
avg_multi[4, ] <- colMeans(table_mononuc_norm2[data_exp > 1 & data_exp <= 3, ])
avg_multi[5, ] <- colMeans(table_mononuc_norm2[data_exp > 3, ])

position <- 1:150

plot(position, avg_multi[1, ], type = "o", col = "gray", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Read count per million mapped read", main = paste0("clone", k))
par(new = T)
plot(position, avg_multi[2, ], type = "o", col = "blue", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Read count per million mapped read", main = paste0("clone", k))
par(new = T)
plot(position, avg_multi[3, ], type = "o", col = "green", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Read count per million mapped read", main = paste0("clone", k))
par(new = T)
plot(position, avg_multi[4, ], type = "o", col = "orange", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Read count per million mapped read", main = paste0("clone", k))
par(new = T)
plot(position, avg_multi[5, ], type = "o", col = "red", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Read count per million mapped read", main = paste0("clone", k))

## save the normalized result
write.table(table_mononuc_norm2, args[11], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)



dev.off()
