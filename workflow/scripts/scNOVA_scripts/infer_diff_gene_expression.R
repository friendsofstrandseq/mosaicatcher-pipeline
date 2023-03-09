args <- commandArgs(trailingOnly = TRUE)

## ---------------------------------------------------------------------------------
## DESeq after filtering out NE from deepCNN
## ---------------------------------------------------------------------------------
library(matrixStats)
library(DESeq2)
library(Rtsne)
library(umap)
library(pheatmap)
library(gplots)

filename = args[8]


prefix <- strsplit(filename, "scNOVA_result_plots")[[1]][1]
prefix <- substring(prefix, 1, nchar(prefix) - 1)
if (nchar(prefix) == 0) {
    prefix <- "."
}
print(prefix)

pdf(filename, width = 11, height = 10)


## 1) Load count matrix

## Load Strand-seq data
data1_GB <- read.table(args[1], sep = "\t", header = T, comment.char = "")

cols <- colnames(data1_GB)
cols <- sub(".sort.mdup.bam", "", cols)
cols <- sub("Deeptool_Genes_for_CNN_", "", cols)
print(cols)
colnames(data1_GB) <- cols

data1_GB_new <- data1_GB[, 4:ncol(data1_GB)]

GB_count <- cbind(data1_GB_new)

## To cluster selected cells in BCLL01
class_label <- rbind(matrix("GM20509", ncol(data1_GB_new), 1))
class_label <- as.matrix(class_label)
class_label_strict_subclone <- read.table(args[2], header = T, sep = "\t")
class_label_sce <- cbind(class_label, as.character(class_label_strict_subclone[, 2]))

cond_label <- class_label_strict_subclone[, 2]

## Sort the order of single-cells
GB_count_name <- as.data.frame(as.matrix(colnames(GB_count)))
GB_count_name$index <- 0
for (j in 1:nrow(GB_count_name)) {
    GB_count_name[j, 1] <- strsplit(GB_count_name[j, 1], ".bam")[[1]][1]
    # GB_count_name[j, 1] <- strsplit(GB_count_name[j, 1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
    GB_count_name[j, 2] <- which(class_label_strict_subclone[, 1] == GB_count_name[j, 1])
}
GB_count <- GB_count[, order(GB_count_name[, 2])]



## Load Gene information
TSS_matrix <- read.table(args[3], header = TRUE, sep = "\t")
GB_matrix <- read.table(args[4], header = TRUE, sep = "\t")


## Sort deeptool result to match with TSS matrix (GB_matrix is same order with deeptool result of GB)
GB_count_sort <- matrix(0, nrow(GB_count), ncol(GB_count))
GB_matrix_sort <- matrix(0, nrow(GB_matrix), ncol(GB_matrix))
colnames(GB_count_sort) <- colnames(GB_count)
colnames(GB_matrix_sort) <- colnames(GB_matrix)
for (i in 1:nrow(GB_count_sort)) {
    GB_count_sort[i, ] <- as.matrix(GB_count[GB_matrix[, 5] == TSS_matrix[i, 5], ])
    GB_matrix_sort[i, ] <- as.matrix(GB_matrix[GB_matrix[, 5] == TSS_matrix[i, 5], ])
    # cat(paste0(i, ' '))
}
GB_matrix_sort <- as.data.frame(GB_matrix_sort)
GB_matrix_sort[, 3] <- as.numeric(as.character(GB_matrix_sort[, 3]))
GB_matrix_sort[, 4] <- as.numeric(as.character(GB_matrix_sort[, 4]))
GB_matrix_sort[, 11] <- as.numeric(as.character(GB_matrix_sort[, 11]))



GB_count_sort_woMT <- GB_count_sort[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]
GB_matrix_sort_woMT <- GB_matrix_sort[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]


## 2) Load CNN prediction and filter out NEs


CNN_result1 <- read.table(args[5], sep = "\t", header = T, comment.char = "")
CNN_result2 <- read.table(args[6], sep = "\t", header = T, comment.char = "")
CNN_result_all <- cbind(CNN_result1[, 8], CNN_result2[, 8])
CNN_result_all_sort <- matrix(0, nrow(TSS_matrix), 2)
for (i in 1:nrow(TSS_matrix)) {
    if (sum(CNN_result1[, 4] == as.character(TSS_matrix[i, 5])) > 0) {
        CNN_result_all_sort[i, ] <- as.matrix(CNN_result_all[CNN_result1[, 4] == as.character(TSS_matrix[i, 5]), ])
        # cat(paste0(i, ' '))
    }
}
colnames(CNN_result_all_sort) <- c("clone1", "clone2")
# filename <- paste0("result_CNN/Expressed_train80_final_result.txt")
filename <- args[9]
write.table(CNN_result_all_sort, file = filename, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

Expressed_pred <- CNN_result_all_sort
Expressed_pred_median <- rowMedians(as.matrix(Expressed_pred[, 1:2]))
Expressed_pred_woMT <- Expressed_pred[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]
Expressed_pred_median_woMT <- Expressed_pred_median[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM"]



CM <- GB_count_sort_woMT
CM_row <- GB_matrix_sort_woMT$name
CM <- round(CM)
CM <- as.data.frame(CM)

rownames(CM) <- as.matrix(CM_row)
for (i in 1:ncol(CM)) {
    CM[, i] <- as.integer(CM[, i])
}




conds <- as.matrix(class_label_strict_subclone[, 2])
conds <- as.factor(conds)
names(conds) <- colnames(CM)



### plot read numbers per sample:
barplot(colSums(CM), col = "sienna", las = 2, cex.names = 0.6, cex.axis = 0.6, ylab = "library size")
mtext(side = 2, line = 4.5, text = "Number of transcript-counted reads", font = 2)


## Filtering out NE from deepCNN
k <- 3

threshold <- (k - 1) * 0.05 ## 0~0.95 (20 levels)
CM_filtered <- CM[is.na(Expressed_pred_woMT[, 1]) == 0 & (Expressed_pred_woMT[, 1] >= threshold | Expressed_pred_woMT[, 2] >= threshold), ]





## 3) Normalization and DESeq

## Make condition column for the sample information
cts <- CM_filtered
coldata <- as.data.frame(conds)
rownames(coldata) <- colnames(cts)
colnames(coldata) <- "condition"


## To extract normalized read count

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
dds <- DESeq(dds)
normcount <- counts(dds, normalized = TRUE)
# write.table(normcount, file = "DESeq_normcount_Strand_GB_LCL_K562_RPE1_2019_19770.txt", sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)


### correlation:
CM.cor <- cor(normcount, method = "spearman")



## general correlation between samples:
CM.dists <- 1 - as.dist(cor(normcount, method = "spearman"))
CM.hc <- hclust(CM.dists)
plot(CM.hc, main = NA, xlab = "Spearman correlation distance")



## Plot PCA

normcount_var <- matrix(0, nrow(normcount), 2)
for (i in 1:nrow(normcount_var)) {
    normcount_var[i, 1] <- sd(normcount[i, ])
}
normcount_var[, 2] <- rank(-1 * normcount_var[, 1])

library("factoextra")
# data_sub_t<-as.data.frame(t(normcount[normcount_var[,1]!=0,])) ##Column: genes, ##Rows: samples
data_sub_t <- as.data.frame(t(normcount[normcount_var[, 2] <= 500, ]))
# data_sub_t<-as.data.frame(t(normcount[normcount_var[,1]!=0 & GB_matrix_sort_woMT[,2]!="chr3",]))
data_pca_scale <- prcomp(data_sub_t, center = TRUE, scale = TRUE)
ind.coord <- data_pca_scale$x

cond_label2 <- as.matrix(cond_label)
data_lab_mat_sub <- cond_label2
data_lab_mat_sub[data_lab_mat_sub == "clone1"] <- "yellowgreen"
data_lab_mat_sub[data_lab_mat_sub == "clone2"] <- "magenta"
plot(ind.coord[, 1], ind.coord[, 2], col = data_lab_mat_sub, pch = 16, xlab = "PC1", ylab = "PC2", cex = 1)
plot(ind.coord[, 1], ind.coord[, 3], col = data_lab_mat_sub, pch = 16, xlab = "PC1", ylab = "PC3", cex = 1)
plot(ind.coord[, 2], ind.coord[, 3], col = data_lab_mat_sub, pch = 16, xlab = "PC2", ylab = "PC3", cex = 1)
plot(ind.coord[, 1], ind.coord[, 4], col = data_lab_mat_sub, pch = 16, xlab = "PC1", ylab = "PC4", cex = 1)


# Eigenvalues
eig <- (data_pca_scale$sdev)^2
# Variances in percentage
variance <- eig * 100 / sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(
    eig = eig, variance = variance,
    cumvariance = cumvar
)
head(eig.decathlon2.active)


# barplot(eig.decathlon2.active[, 2], names.arg=1:nrow(eig.decathlon2.active),
barplot(eig.decathlon2.active[1:10, 2],
    names.arg = 1:10,
    main = "Variances",
    xlab = "Principal Components",
    ylab = "Percentage of variances",
    col = "steelblue"
)
# Add connected line segments to the plot
lines(
    x = 1:10,
    eig.decathlon2.active[1:10, 2],
    type = "b", pch = 19, col = "red"
)



d <- dist(ind.coord[, 1:10])
# print(ind.coord)
# print(nrow(ind.coord))
# print(floor((nrow(ind.coord) - 1) / 3))
# print(d)
# print(nrow(d))
# print(floor((nrow(d) - 1) / 3))
set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random
# tsne_out <- Rtsne(d, is_distance = TRUE, perplexity = 10, verbose = TRUE)
tsne_out <- Rtsne(d, is_distance = TRUE, perplexity = floor((nrow(ind.coord) - 1) / 3), verbose = TRUE)
plot(tsne_out$Y, pch = 16, xlab = "t-SNE1", ylab = "t-SNE2", cex = 1, col = data_lab_mat_sub)
# write.table(tsne_out$Y, file = "/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_ATAC/Active_X_haplo_analysis/LCL_GM20509/output_tsne_500genes.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


umap_out <- umap(ind.coord[, 1:10])
plot(umap_out$layout, pch = 16, xlab = "UMAP1", ylab = "UMAP2", cex = 1, col = data_lab_mat_sub)
# write.table(umap_out$layout, file = "/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_ATAC/Active_X_haplo_analysis/LCL_GM20509/output_umap_500genes.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)



## To identify DEGs (Wald test, by default)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
dds <- DESeq(dds)
normcount <- counts(dds, normalized = TRUE)
res <- results(dds, contrast = c("condition", "clone2", "clone1"))


## Sort DESeq result to match and annotate TSS_matrix
res_sort <- as.data.frame(matrix(1, nrow(GB_matrix_sort), ncol(res))) # If there's no datapoint, FDR will be 1 by default
for (i in 1:nrow(GB_matrix_sort)) {
    if (sum(rownames(res) == GB_matrix_sort[i, 5]) > 0) {
        res_sort[i, ] <- as.matrix(res[rownames(res) == GB_matrix_sort[i, 5], ])
    }
    # cat(paste0(i, ' '))
}
colnames(res_sort) <- colnames(res)


res_sort_thres3 <- res_sort





## 4) Single-cell heatmap generation

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
res_sort_woMT <- res_sort[TSS_matrix[, 2] != "chrY" & TSS_matrix[, 2] != "chrM", ]


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

if (sum(input_matrix_sort_woMT$blacklist == 0 & res_sort_woMT$padj < 0.1 & is.na(res_sort_woMT$padj) == 0) > 1) {
    breaksList <- seq(-1.5, 1.5, by = 0.1)
    breaksList <- append(breaksList, 2)
    breaksList <- append(breaksList, -2, 0)
    mycol <- colorpanel(n = length(breaksList) - 1, low = "blue", mid = "white", high = "red")
    res <- pheatmap(normlogt[class_label == "clone1" | class_label == "clone2", input_matrix_sort_woMT$blacklist == 0 & res_sort_woMT$padj < 0.1 & is.na(res_sort_woMT$padj) == 0], show_rownames = F, show_colnames = T, cluster_cols = T, cluster_rows = T, scale = "column", col = mycol, breaks = breaksList, clustering_distance_rows = "euclidean", cex = 0.8, annotation_row = row_annotation, annotation_colors = anno_colors, clustering_method = "ward.D")


    ## Without clustering row
    normlogt_sort <- rbind(normlogt[class_label == "clone1", ], normlogt[class_label == "clone2", ])
    row_annotation_sort <- rbind(as.matrix(row_annotation[class_label == "clone1", ]), as.matrix(row_annotation[class_label == "clone2", ]))
    rownames(row_annotation_sort) <- rownames(normlogt_sort)
    colnames(row_annotation_sort) <- "Var1"
    pheatmap(normlogt_sort[, input_matrix_sort_woMT$blacklist == 0 & res_sort_woMT$padj < 0.1 & is.na(res_sort_woMT$padj) == 0], show_rownames = F, show_colnames = T, cluster_cols = T, cluster_rows = F, scale = "column", col = mycol, breaks = breaksList, clustering_distance_rows = "euclidean", cex = 0.8, annotation_row = row_annotation, annotation_colors = anno_colors, clustering_method = "ward.D")
}

res_sort_hit <- cbind(res_sort, input_matrix$blacklist, Expressed_pred_median)
res_sort_hit$Hit <- 0
res_sort_hit[input_matrix$blacklist == 0 & res_sort$padj < 0.1 & is.na(res_sort$padj) == 0, ncol(res_sort_hit)] <- 1

write.table(cbind(GB_matrix_sort[, 2:5], res_sort_hit), file = paste0(prefix, "/scNOVA_result/Result_scNOVA_infer_expression_table.txt"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

dev.off()
