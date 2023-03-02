args <- commandArgs(trailingOnly = TRUE)

## ---------------------------------------------------------------------------------
## DESeq after filtering out NE from deepCNN (final checking!! 20210723)
## ---------------------------------------------------------------------------------


## ----------------------------------------
## 1) Load Strand-seq count matrix
## ----------------------------------------

library(matrixStats)
library(DESeq2)
library(Rtsne)
library(umap)
library(pheatmap)
library(gplots)
library(fitdistrplus)
library(doParallel)
library(foreach)

pdf(args[9], width = 11, height = 10)

## 1) Load count matrix


## Load Strand-seq data (Condition 1)
data1_GB <- read.table(args[1], sep = "\t", header = T, comment.char = "")

cols <- colnames(data1_GB)
cols <- sub(".sort.mdup.bam", "", cols)
cols <- sub("Deeptool_Genes_for_CNN_", "", cols)
print(cols)
colnames(data1_GB) <- cols

data1_GB_new <- data1_GB[, 4:ncol(data1_GB)]


print("A")


GB_count <- cbind(data1_GB_new)
# GB_count <- GB_count[data1_GB[,1]!="chrY" & data1_GB[,1]!="chrM",]##ChrX will be included in the analysis!!

## To cluster selected cells in BCLL01
class_label <- rbind(matrix("DEA5", ncol(data1_GB_new), 1))
class_label <- as.matrix(class_label)
class_label_strict_subclone <- read.table(args[2], header = T, sep = "\t")
class_label_sce <- cbind(class_label, as.character(class_label_strict_subclone[, 2]))

cond_label <- class_label_strict_subclone[, 2]


print("B")

## Sort the order of single-cells
GB_count_name <- as.data.frame(as.matrix(colnames(GB_count)))
GB_count_name$index <- 0
for (j in 1:nrow(GB_count_name)) {
    GB_count_name[j, 1] <- strsplit(GB_count_name[j, 1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
    GB_count_name[j, 2] <- which(class_label_strict_subclone[, 1] == GB_count_name[j, 1])
}
GB_count <- GB_count[, order(GB_count_name[, 2])]

print("C")

## Load Gene information

TSS_matrix <- read.table(args[3], header = TRUE, sep = "\t")
GB_matrix <- read.table(args[4], header = TRUE, sep = "\t")

print("C2")


# GB_matrix <- GB_matrix[data1_GB[,1]!="chrY" & data1_GB[,1]!="chrM",]

## Sort deeptool result to match with TSS matrix (GB_matrix is same order with deeptool result of GB)
GB_count_sort <- matrix(0, nrow(GB_count), ncol(GB_count))
# GB_count_sort <- matrix(0, 100, ncol(GB_count))
GB_matrix_sort <- matrix(0, nrow(GB_matrix), ncol(GB_matrix))
# GB_matrix_sort <- matrix(0, 100, ncol(GB_matrix))
colnames(GB_count_sort) <- colnames(GB_count)
colnames(GB_matrix_sort) <- colnames(GB_matrix)
print("C3")

# pb <- progress_bar$new(total = nrow(GB_count_sort))
print(nrow(GB_count_sort))
# for (i in 1:nrow(GB_count_sort)) {
# pb$tick()
# GB_count_sort[i, ] <- as.matrix(GB_count[GB_matrix[, 5] == TSS_matrix[i, 5], ])
GB_count_sort <- GB_count[match(TSS_matrix[, 5], GB_matrix[, 5]), ]
GB_matrix_sort <- GB_matrix[match(TSS_matrix[, 5], GB_matrix[, 5]), ]

# print(head(GB_count_sort))
# print(nrow(GB_count_sort))
# GB_matrix_sort[i, ] <- as.matrix(GB_matrix[GB_matrix[, 5] == TSS_matrix[i, 5], ])
# cat(paste0(i, ' '))
# }
print("C4")


GB_matrix_sort <- as.data.frame(GB_matrix_sort)
GB_matrix_sort[, 3] <- as.numeric(as.character(GB_matrix_sort[, 3]))
GB_matrix_sort[, 4] <- as.numeric(as.character(GB_matrix_sort[, 4]))
GB_matrix_sort[, 11] <- as.numeric(as.character(GB_matrix_sort[, 11]))
# write.table(GB_count_sort, "/Users/jeong/Desktop/GB_count_sort.txt", row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)

print(head(GB_count_sort))
print(head(GB_matrix_sort))
# stop()
print("D")


## To annotate cells with and without mSV
conds <- as.matrix(class_label_strict_subclone[, 2])
conds <- as.factor(conds)
names(conds) <- colnames(GB_count_sort)


cond_label <- conds
Y <- as.matrix(cond_label)
Y[Y[, 1] == "clone2", 1] <- 1
Y[Y[, 1] == "clone1", 1] <- 0
Y <- as.matrix(as.numeric(Y))
Y <- as.data.frame(Y)
colnames(Y) <- "SV"
Y$none <- 0
Y[, 2] <- 1 - Y[, 1]
# write.table(Y, "/Users/jeong/Desktop/Y.txt", row.names = TRUE, col.names = TRUE, sep="\t", quote = FALSE)

print("E")

## ----------------------------------------
## 2) RPM normalization, Load CNN prediction and filter out NEs
## ----------------------------------------

## Load CNN prediction and filter out NEs



prefix <- strsplit(output_filename, "scNOVA_result")[[1]][1]
prefix <- substring(prefix, 1, nchar(prefix) - 1)
if (nchar(prefix) == 0) {
    prefix <- "."
}
print(prefix)

Expressed_pred <- read.table(args[10]), sep = "\t", header = T, comment.char = "")
Expressed_pred_median <- rowMedians(as.matrix(Expressed_pred[, 1:2]))

print("F")


## Filtering out NE from deepCNN
k <- 3
threshold <- (k - 1) * 0.05 ## 0~0.95 (20 levels)
CNN_filter <- (is.na(Expressed_pred[, 1]) == 0 & (Expressed_pred[, 1] >= threshold | Expressed_pred[, 2] >= threshold))

print("G")


data_GB_RPM <- matrix(0, nrow(GB_count_sort), ncol(GB_count_sort))

for (i in 1:ncol(GB_count_sort)) {
    data_GB_RPM[, i] <- GB_count_sort[, i] * 1000000 / sum(GB_count_sort[, i])
}
print("H")


hist(data_GB_RPM[, 1], 200)


data_GB_RPM_log <- log2(data_GB_RPM + 1)
data_GB_RPM_log_CV <- matrix(0, nrow(GB_count_sort), 1)

for (i in 1:nrow(data_GB_RPM_log)) {
    data_GB_RPM_log_CV[i, 1] <- sd(data_GB_RPM_log[i, ]) / mean(data_GB_RPM_log[i, ])
}

print("I")


data_GB_RPM_log_sub <- data_GB_RPM_log[is.na(data_GB_RPM_log_CV) == 0 & data_GB_RPM_log_CV > 0 & CNN_filter > 0, ]
X_GB_RPM_log <- t(data_GB_RPM_log_sub)
# write.table(data_GB_RPM_log_CV, file = "/Users/jeong/Desktop/data_GB_RPM_log_CV.txt", sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

print("J")


## ----------------------------------------
## 3) Make X and Y matrix and model building
## ----------------------------------------


Y <- as.matrix(Y)


## ----------------------------------------
## 4) Autoscaling X and Y matrix
## ----------------------------------------

source("workflow/scripts/scNOVA_scripts/script_PLSDA/auto_R.R")
X_auto <- auto_R(X_GB_RPM_log)
Y_auto <- auto_R(Y)


print("K")


## ----------------------------------------
## 5) Run PLS-DA
## ----------------------------------------


# First pls model using all features
if (length(conds) > 20) {
    lv <- 20
}
if (length(conds) <= 20) {
    lv <- (length(conds) - 1)
}

source("workflow/scripts/scNOVA_scripts/script_PLSDA/pls_R_scNOVA.R")
result_pls_all <- pls_R(X_auto, Y_auto, lv)
# Plot LV1 and LV2
# data_lab_mat_sub <- as.matrix(conds)
# data_lab_mat_sub[data_lab_mat_sub == "clone1"] <- "yellowgreen"
# data_lab_mat_sub[data_lab_mat_sub == "clone2"] <- "magenta"
# data_lab_mat_sub[data_lab_mat_sub == "clone3"] <- "blue"
# plot(result_pls_all$pls_t[, 1], result_pls_all$pls_t[, 2], col = data_lab_mat_sub, pch = 16, xlab = "LV1", ylab = "LV2", cex = 1)

print("L")



# cl <- makeCluster(20) # create a cluster with all available cores except one
# registerDoParallel(cl)

# iterations <- lv
# pb <- txtProgressBar(max = iterations, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# evaluation <- foreach(aj = 2:lv, .combine = rbind, .options.snow = opts) %dopar% {
#     tmp <- matrix(0, nrow(X_auto), 1)
#     for (ai in 1:nrow(X_auto)) {
#         X_auto_train <- X_auto[-ai, ]
#         Y_auto_train <- Y_auto[-ai, ]
#         X_auto_test <- X_auto[ai, ]
#         Y_auto_test <- Y_auto[ai, ]
#         source("workflow/scripts/scNOVA_scripts/script_PLSDA/Pred_PLS_R_scNOVA.R")
#         result_ypred <- Pred_PLS_R(X_auto_train, Y_auto_train, X_auto_test, aj)
#         tmp[ai, 1] <- sum(which.max(Y_auto_test) == which.max(result_ypred))
#     }
#     tmp
# }
# stopCluster(cl)

cl <- makeCluster(64, outfile = "cluster.log")
registerDoParallel(cl)


# Pre-allocate evaluation matrix
evaluation <- matrix(0, nrow(X_auto), lv)

# Pre-allocate train and test sets
X_auto_train <- matrix(0, nrow(X_auto) - 1, ncol(X_auto))
Y_auto_train <- matrix(0, nrow(Y_auto) - 1, ncol(Y_auto))
X_auto_test <- matrix(0, 1, ncol(X_auto))
Y_auto_test <- matrix(0, 1, ncol(Y_auto))

source("workflow/scripts/scNOVA_scripts/script_PLSDA/Pred_PLS_R_scNOVA.R")


# Set progress bar options
opts <- list(progress = "text", verbose = FALSE)

for (aj in 2:lv) {
    print(aj)
    # Prepare train and test sets
    X_auto_train <- X_auto[-1, ]
    Y_auto_train <- Y_auto[-1, ]
    X_auto_test <- X_auto[1, ]
    Y_auto_test <- Y_auto[1, ]

    # Define foreach loop
    results <- foreach(ai = 1:nrow(X_auto), .combine = rbind, .options.snow = opts) %dopar% {
        result_ypred <- Pred_PLS_R(X_auto_train[-ai, ], Y_auto_train[-ai, ], X_auto_test, aj)
        sum(which.max(Y_auto_test) == which.max(result_ypred))
    }

    # Store results in evaluation matrix
    evaluation[, aj] <- results
}

# Stop the parallel backend
stopCluster(cl)


# # Leave-One-Out cross validation using all features
# evaluation <- matrix(0, nrow(X_auto), lv)
# for (aj in 2:lv) {
#     for (ai in 1:nrow(X_auto)) {
#         X_auto_train <- X_auto[-ai, ]
#         Y_auto_train <- Y_auto[-ai, ]
#         X_auto_test <- X_auto[ai, ]
#         Y_auto_test <- Y_auto[ai, ]
#         source("workflow/scripts/scNOVA_scripts/script_PLSDA/Pred_PLS_R_scNOVA.R")
#         result_ypred <- Pred_PLS_R(X_auto_train, Y_auto_train, X_auto_test, aj)
#         evaluation[ai, aj] <- sum(which.max(Y_auto_test) == which.max(result_ypred))
#     }
#     cat(paste0(aj, " "))
# }
plot(colMeans(evaluation), type = "l", xlab = "num of LV", ylab = "Accuracy", ylim = c(0, 1))


print("M")



# Feature selection using VIP
lv_for_vip <- which.max(colMeans(evaluation))
w <- result_pls_all$pls_w[, 1:lv_for_vip]
r2 <- result_pls_all$pls_ssq[1:lv_for_vip, 4]
source("workflow/scripts/scNOVA_scripts/script_PLSDA/vip_R_v2_scNOVA.R")
result_vip_all <- vip_R_v2(w, r2)


# Generation of null distribution of VIP
source("workflow/scripts/scNOVA_scripts/script_PLSDA/vip_null_R.R")
perm <- 100
result_vip_all_null <- vip_null_R(X_auto, Y_auto, perm, lv_for_vip)
hist(result_vip_all_null, 200, xlab = "VIP", ylab = "Frequency", main = "empirical distribution of VIP")


print("N")

## ----------------------------------------
## 6) Calculate p-value and FDR
## ----------------------------------------



data_null <- result_vip_all_null
data_vip <- result_vip_all

data_null <- as.matrix(data_null)
data_vip <- as.matrix(data_vip)


fno <- fitdist(as.vector(data_null), "norm")
denscomp(fno)
cdfcomp(fno)
summary(fno)


print("O")


fno_pval <- 1 - pnorm(data_vip, mean = fno$estimate[1], sd = fno$estimate[2])
fno_fdr <- p.adjust(fno_pval, method = "BH")
fno_fdr2 <- p.adjust(fno_pval, method = "bonferroni")

result_table <- cbind(fno_pval, fno_fdr, fno_fdr2)
colnames(result_table) <- c("fno_pval", "fno_fdr", "fno_fdr2")
# write.table(result_table, file = "/Users/jeong/Desktop/result_table_CV0_CNNfilter_DEA5.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


GB_matrix_sort_test <- GB_matrix_sort[is.na(data_GB_RPM_log_CV) == 0 & data_GB_RPM_log_CV > 0 & CNN_filter > 0, ]

result_table_all <- matrix(0, nrow(GB_matrix_sort), 4)
for (i in 1:nrow(GB_matrix_sort_test)) {
    if (sum(GB_matrix_sort[, 5] == GB_matrix_sort_test[i, 5]) > 0) {
        result_table_all[which(GB_matrix_sort[, 5] == GB_matrix_sort_test[i, 5]), 1] <- result_table[i, 1]
        result_table_all[which(GB_matrix_sort[, 5] == GB_matrix_sort_test[i, 5]), 2] <- result_table[i, 2]
        result_table_all[which(GB_matrix_sort[, 5] == GB_matrix_sort_test[i, 5]), 3] <- result_table[i, 3]
        result_table_all[which(GB_matrix_sort[, 5] == GB_matrix_sort_test[i, 5]), 4] <- 1
    }
}
colnames(result_table_all) <- c("fno_pval", "fno_fdr", "fno_fdr2", "test")


data_GB_RPM_log_FC <- rowMeans(data_GB_RPM_log[, Y[, 1] == 1]) - rowMeans(data_GB_RPM_log[, Y[, 1] == 0])
result_table_all_annot <- cbind(GB_matrix_sort[, 1:5], result_table_all)
result_table_all_annot$log2FC <- as.matrix(data_GB_RPM_log_FC)
result_table_all_annot$Expressed <- CNN_filter


SV_affected <- read.table(args[7], header = TRUE, sep = "\t")
result_table_all_annot$blacklist <- SV_affected$blacklist
result_table_all_annot$Hit <- 0
result_table_all_annot[result_table_all_annot$test == 1 & result_table_all_annot$fno_fdr < 0.1 & result_table_all_annot$blacklist == 0, ncol(result_table_all_annot)] <- 1


write.table(result_table_all_annot, file = args[8], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)



## ----------------------------------------
## 7) Heatmap of selected genes
## ----------------------------------------


GB_count_sort_woMT <- GB_count_sort[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]
GB_matrix_sort_woMT <- GB_matrix_sort[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]

CM <- GB_count_sort_woMT
CM_row <- GB_matrix_sort_woMT$name
CM <- round(CM)
CM <- as.data.frame(CM)

rownames(CM) <- as.matrix(CM_row)
for (i in 1:ncol(CM)) {
    CM[, i] <- as.integer(CM[, i])
}

class_label_strict_subclone_v2 <- as.data.frame(as.matrix(conds))




## Make condition column for the sample information
cts <- CM
coldata <- as.data.frame(conds)
rownames(coldata) <- colnames(cts)
colnames(coldata) <- "condition"




## To extract normalized read count
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
dds <- DESeq(dds)
normcount <- counts(dds, normalized = TRUE)


input_matrix <- read.table(args[7], sep = "\t", header = T, comment.char = "")
input_matrix_sort_woMT <- input_matrix[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]
res_sort_woMT <- result_table_all_annot[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]


normlog <- log2(normcount + 1)
normlogt <- t(normlog)

class_label <- as.character(class_label_strict_subclone[, 2])
class_label2 <- as.matrix(class_label)
class_label2[class_label == "clone1"] <- "yellowgreen"
class_label2[class_label == "clone2"] <- "magenta"


Var1 <- c("yellowgreen", "magenta")
names(Var1) <- c("yellowgreen", "magenta")
anno_colors <- list(Var1 = Var1)
row_annotation <- as.data.frame(class_label2)
rownames(row_annotation) <- colnames(normcount)
colnames(row_annotation) <- "Var1"


## Generate heatmap

if (sum(res_sort_woMT$Hit == 1) > 1) {
    breaksList <- seq(-1.5, 1.5, by = 0.1)
    breaksList <- append(breaksList, 2)
    breaksList <- append(breaksList, -2, 0)
    mycol <- colorpanel(n = length(breaksList) - 1, low = "blue", mid = "white", high = "red")
    res <- pheatmap(normlogt[class_label == "clone1" | class_label == "clone2", res_sort_woMT$Hit == 1], show_rownames = F, show_colnames = T, cluster_cols = T, cluster_rows = T, scale = "column", col = mycol, breaks = breaksList, clustering_distance_rows = "euclidean", cex = 0.8, annotation_row = row_annotation, annotation_colors = anno_colors, clustering_method = "ward.D")


    ## Without clustering row
    normlogt_sort <- rbind(normlogt[class_label == "clone1", ], normlogt[class_label == "clone2", ])
    row_annotation_sort <- rbind(as.matrix(row_annotation[class_label == "clone1", ]), as.matrix(row_annotation[class_label == "clone2", ]))
    rownames(row_annotation_sort) <- rownames(normlogt_sort)
    colnames(row_annotation_sort) <- "Var1"
    pheatmap(normlogt_sort[, res_sort_woMT$Hit == 1], show_rownames = F, show_colnames = T, cluster_cols = T, cluster_rows = F, scale = "column", col = mycol, breaks = breaksList, clustering_distance_rows = "euclidean", cex = 0.8, annotation_row = row_annotation, annotation_colors = anno_colors, clustering_method = "ward.D")
}




dev.off()
