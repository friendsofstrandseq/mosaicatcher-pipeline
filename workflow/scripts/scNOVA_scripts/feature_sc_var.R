args <- commandArgs(trailingOnly = TRUE)

output_filename <- args[8]

pdf(output_filename, width = 11, height = 10)



prefix <- strsplit(output_filename, "scNOVA_result")[[1]][1]
prefix <- substring(prefix, 1, nchar(prefix) - 1)
if (nchar(prefix) == 0) {
    prefix <- "."
}
print(prefix)


subclone_list <- read.table(args[1], header = T, sep = "\t", comment.char = "")
subclone_count <- length(unique(subclone_list[, 2]))
subclone_uniq <- sort(unique(subclone_list[, 2]))

Deeptool_result_final <- read.table(args[2], header = TRUE, sep = "\t", comment.char = "")

cols <- colnames(Deeptool_result_final)
cols <- sub(".sort.mdup.bam", "", cols)
cols <- sub("Deeptool_Genes_for_CNN_", "", cols)
print(cols)
colnames(Deeptool_result_final) <- cols

Deeptool_result_final_new_all <- Deeptool_result_final[, 5:ncol(Deeptool_result_final)]

print(head(Deeptool_result_final_new_all))

print("A")
## Sort the order of single-cells
Deeptool_result_name <- as.data.frame(as.matrix(colnames(Deeptool_result_final_new_all)))
# print("B")

Deeptool_result_name$index <- 0
for (j in 1:nrow(Deeptool_result_name)) {
    print(j)
    # print("C")
    Deeptool_result_name[j, 1] <- strsplit(Deeptool_result_name[j, 1], ".bam")[[1]][1]
    # Deeptool_result_name[j, 1] <- strsplit(Deeptool_result_name[j, 1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
    # print("D")
    Deeptool_result_name[j, 2] <- which(subclone_list[, 1] == Deeptool_result_name[j, 1])
}
Deeptool_result_final_new_all <- Deeptool_result_final_new_all[, order(Deeptool_result_name[, 2])]
print("B")


# for (k in 1:subclone_count){
clone_name <- strsplit(args[7], paste0(prefix, "/scNOVA_result/Features_reshape_"))[[1]][2]
clone_name <- strsplit(clone_name, "_orientation_CN_correct0.txt")[[1]][1]
filename <- paste0(prefix, "/scNOVA_result/Features_reshape_", clone_name, "_Resid_orientation.txt")
filename_IQR <- paste0(prefix, "/scNOVA_result/Features_reshape_", clone_name, "_Resid_orientation_IQR.txt")
Deeptool_result_final_new <- Deeptool_result_final_new_all[, (subclone_list[, 2] == clone_name)]
Ref_bed_annot <- read.table(args[3], header = T, sep = " ", comment.char = "")
Ref_bed_annot <- Ref_bed_annot[Ref_bed_annot[, 1] != "chrM", ]

print("C")


# Deeptool_result_final_metric <- matrix(0, nrow(Deeptool_result_final), 3)
# for (i in 1:nrow(Deeptool_result_final_new)) {
#     Deeptool_result_final_metric[i, 1] <- mean(as.numeric(Deeptool_result_final_new[i, ])) ## mean
#     Deeptool_result_final_metric[i, 2] <- sd(as.numeric(Deeptool_result_final_new[i, ])) ## standard deviation
#     Deeptool_result_final_metric[i, 3] <- sd(as.numeric(Deeptool_result_final_new[i, ])) / mean(as.numeric(Deeptool_result_final_new[i, ])) ## coefficient of variation
#     # cat(paste0(i, ' '))
# }
Deeptool_result_final_metric <- t(rbind(apply(Deeptool_result_final_new, 1, mean), apply(Deeptool_result_final_new, 1, sd), apply(Deeptool_result_final_new, 1, sd) / apply(Deeptool_result_final_new, 1, mean)))


print("D")

TSS_matrix <- read.table(args[4], header = TRUE, sep = "\t")
# Deeptool_result_final_metric <- Deeptool_result_final_metric
plot(log10(Deeptool_result_final_metric[, 1] + 1), Deeptool_result_final_metric[, 3]^2, xlab = "log10(mean+1)", ylab = "CV^2")

Deeptool_result_final_metric_sub <- Deeptool_result_final_metric[is.na(Deeptool_result_final_metric[, 3]) == 0, ]
Deeptool_result_final_metric_sub_lab <- which(is.na(Deeptool_result_final_metric[, 3]) == 0)

print("E")


if (IQR(log10(Deeptool_result_final_metric_sub[, 1] + 1)) != 0) {
    fit1 <- smooth.spline(log10(Deeptool_result_final_metric_sub[, 1] + 1), Deeptool_result_final_metric_sub[, 3]^2, df = 16)
}
if (IQR(log10(Deeptool_result_final_metric_sub[, 1] + 1)) == 0) {
    fit1 <- smooth.spline(log10(Deeptool_result_final_metric_sub[, 1] + 1), Deeptool_result_final_metric_sub[, 3]^2, df = 16, tol = 1e-6 * 0.146128)
}
write.table(IQR(log10(Deeptool_result_final_metric_sub[, 1] + 1)), filename_IQR, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

print("F")


lines(fit1, col = "red", lwd = 2)
a <- predict(fit1, log10(Deeptool_result_final_metric_sub[, 1] + 1))
Resid <- (Deeptool_result_final_metric_sub[, 3])^2 - a$y
plot(log10(Deeptool_result_final_metric_sub[, 1] + 1), Resid)

result <- matrix(0, nrow(Deeptool_result_final_metric), 1)
for (i in 1:nrow(Deeptool_result_final_metric_sub)) {
    result[Deeptool_result_final_metric_sub_lab[i], 1] <- Resid[i]
}
print("G")


library(pracma)
CNN_features_annot <- read.table(args[5], header = T, sep = "\t", comment.char = "")
table_original <- t(Reshape(result[, 1], 150, (19770 - 13)))
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
write.table(table_final, args[9], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
print("H")

FPKM <- read.table(args[6], header = T, sep = "\t", comment.char = "")
TSS_matrix_woM <- TSS_matrix[TSS_matrix[, 2] != "chrM", ]
FPKM_woM <- FPKM[TSS_matrix[, 2] != "chrM", ]
table_mononuc_var <- table_final
print("I")

plot(c(1:150), colMeans(table_mononuc_var, na.rm = T), type = "o")
data_exp <- rowMeans(FPKM_woM[, 1:9]) ## This parameter need to be changed for each samples
avg_multi <- matrix(0, 5, 150)
avg_multi[1, ] <- colMeans(table_mononuc_var[data_exp == 0, ], na.rm = T)
avg_multi[2, ] <- colMeans(table_mononuc_var[data_exp > 0 & data_exp < 0.1, ], na.rm = T)
avg_multi[3, ] <- colMeans(table_mononuc_var[data_exp > 0.1 & data_exp <= 1, ], na.rm = T)
avg_multi[4, ] <- colMeans(table_mononuc_var[data_exp > 1 & data_exp <= 3, ], na.rm = T)
avg_multi[5, ] <- colMeans(table_mononuc_var[data_exp > 3, ], na.rm = T)
print("J")

position <- 1:150

plot(position, avg_multi[1, ], type = "o", col = "gray", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Resid", main = clone_name)
par(new = T)
plot(position, avg_multi[2, ], type = "o", col = "blue", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Resid", main = clone_name)
par(new = T)
plot(position, avg_multi[3, ], type = "o", col = "green", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Resid", main = clone_name)
par(new = T)
plot(position, avg_multi[4, ], type = "o", col = "orange", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Resid", main = clone_name)
par(new = T)
plot(position, avg_multi[5, ], type = "o", col = "red", cex = 0.1, lwd = 3, ylim = range(avg_multi), xlab = "genebody", ylab = "Resid", main = clone_name)

# }

dev.off()
