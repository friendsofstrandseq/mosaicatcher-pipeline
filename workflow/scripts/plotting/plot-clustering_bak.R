
plot.clustering <- function(inputfile, bin.bed.filename, position.outputfile, chromosome.outputfile) {
    genome_bins <- read.table(bin.bed.filename, sep = "\t", header = F, comment.char = "")
    list_directory <- dir("./", full.names = TRUE)

    pdf(position.outputfile, width = 11, height = 10)


    data1 <- read.table(inputfile, sep = "\t", header = T, comment.char = "")

    data1$color <- 0
    ash12rainbow <- c("#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411")
    sv_call_name <- c("del_h1", "del_h2", "del_hom", "dup_h1", "dup_h2", "dup_hom", "inv_h1", "inv_h2", "inv_hom", "idup_h1", "idup_h2", "complex")

    for (j in 1:length(sv_call_name)) {
        tmp <- which(data1[, 9] == sv_call_name[j])
        data1[tmp, 15] <- j
    }

    data1_pos <- data1[, 1:3]
    data1_pos_uniq <- unique(data1_pos)
    data1_pos_uniq_sort <- data1_pos_uniq[data1_pos_uniq$chrom == "chr1", ]
    chrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

    for (i in 2:length(chrom)) {
        data1_pos_uniq_sort <- rbind(data1_pos_uniq_sort, data1_pos_uniq[data1_pos_uniq$chrom == chrom[i], ])
    }

    data1_pos_uniq_sort$posind <- c(1:nrow(data1_pos_uniq_sort))


    data1_cell <- data1[, 5]
    data1_cell_uniq <- unique(data1_cell)
    data1_cell_uniq_sort <- as.matrix(sort(data1_cell_uniq))

    result <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
    result_sv <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))

    for (i in 1:nrow(data1)) {
        pos_ind <- which(data1_pos_uniq_sort[, 1] == data1[i, 1] & data1_pos_uniq_sort[, 2] == data1[i, 2] & data1_pos_uniq_sort[, 3] == data1[i, 3])
        cell_ind <- which(data1_cell_uniq_sort[, 1] == data1[i, 5])
        result[cell_ind, pos_ind] <- data1[i, 13]
        result_sv[cell_ind, pos_ind] <- data1[i, 15]
    }

    rownames(result) <- data1_cell_uniq_sort
    colnames(result) <- data1_pos_uniq_sort$posind
    result[result == Inf] <- max(c(1, result[result != Inf]))
    rownames(result_sv) <- data1_cell_uniq_sort
    colnames(result_sv) <- data1_pos_uniq_sort$posind

    ## Add chromosome color code to heatmap
    data1_pos_uniq_sort$color <- "magenta"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chrom == "chr2", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr4", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr6", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr8", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr10", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr12", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr14", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr16", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr18", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr20", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chr22", 5] <- "purple"
    data1_pos_uniq_sort[data1_pos_uniq_sort$chr == "chrY", 5] <- "purple"


    library(pheatmap)
    library(gplots)
    library(ComplexHeatmap)

    breaksList <- seq(4, 30, by = 0.1)
    breaksList <- append(breaksList, max(result))
    breaksList <- append(breaksList, -1, 0)
    breaksList <- unique(sort(breaksList))
    mycol <- colorpanel(n = length(breaksList) - 1, low = "lightgoldenrodyellow", mid = "darkorange", high = "firebrick2")


    Var1 <- c("magenta", "purple")
    names(Var1) <- c("magenta", "purple")
    anno_colors <- list(Var1 = Var1)
    col_annotation <- as.data.frame(data1_pos_uniq_sort$color)
    rownames(col_annotation) <- colnames(result)
    colnames(col_annotation) <- "Var1"


    res <- pheatmap(result, border_color = NA, show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = T, clustering_method = "ward.D", scale = "none", col = mycol, main = inputfile, annotation_col = col_annotation, annotation_colors = anno_colors, breaks = breaksList, annotation_legend = FALSE)
    # res <- pheatmap(result, border_color = NA, show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = T, clustering_method = "ward.D", scale = "none", col = mycol, cex = 0.7, main = inputfile, annotation_col = col_annotation, annotation_colors = anno_colors, breaks = breaksList, annotation_legend = FALSE)

    chr_name <- matrix("chr", ncol(result), 1)
    for (i in 1:ncol(result)) {
        chr_name[i, 1] <- paste0(data1_pos_uniq_sort[i, 1], "_", data1_pos_uniq_sort[i, 2], "_", data1_pos_uniq_sort[i, 3])
    }

    colnames(result_sv) <- chr_name
    colors <- structure(c("white", "#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411"), names = c(0:12))

    sv_list <- c("none", "del_h1", "del_h2", "del_hom", "dup_h1", "dup_h2", "dup_hom", "inv_h1", "inv_h2", "inv_hom", "idup_h1", "idup_h2", "complex")
    sv_tmp <- matrix(0, 13, 1)
    for (i in 1:13) {
        sv_tmp[i, 1] <- sum(result_sv == (i - 1))
    }
    sv_list_sub <- sv_list[sv_tmp > 0]

    mat <- as.data.frame(result_sv[res$tree_row$order, ])
    ha_column <- HeatmapAnnotation(
        df = data.frame(type1 = data1_pos_uniq_sort$color),
        col = list(type1 = c("magenta" = "magenta", "purple" = "purple"))
    )
    ht1 <- Heatmap(mat, name = "", col = colors, heatmap_legend_param = list(labels = sv_list_sub), cluster_rows = FALSE, cluster_columns = FALSE, column_title = inputfile, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 1), column_title_gp = gpar(fontsize = 7, fontface = "bold"), top_annotation = ha_column)
    draw(ht1, show_annotation_legend = FALSE)

    dev.off()

    pdf(chromosome.outputfile, width = 11, height = 5)

    data1 <- read.table(inputfile, sep = "\t", header = T, comment.char = "")

    data1$color <- 0
    ash12rainbow <- c("#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411")
    sv_call_name <- c("del_h1", "del_h2", "del_hom", "dup_h1", "dup_h2", "dup_hom", "inv_h1", "inv_h2", "inv_hom", "idup_h1", "idup_h2", "complex")
    for (j in 1:length(sv_call_name)) {
        tmp <- which(data1[, 9] == sv_call_name[j])
        data1[tmp, 15] <- j
    }
    colors <- structure(c("white", "#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411"), names = c(0:12))


    data1_pos <- data1[, 1:3]
    data1_pos_uniq <- unique(data1_pos)
    data1_pos_uniq_sort <- data1_pos_uniq[data1_pos_uniq$chrom == "chr1", ]
    chrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

    for (i in 2:length(chrom)) {
        data1_pos_uniq_sort <- rbind(data1_pos_uniq_sort, data1_pos_uniq[data1_pos_uniq$chrom == chrom[i], ])
    }
    
    data1_pos_uniq_sort$posind <- c(1:nrow(data1_pos_uniq_sort))

    data1_cell <- data1[, 5]
    data1_cell_uniq <- unique(data1_cell)
    data1_cell_uniq_sort <- as.matrix(sort(data1_cell_uniq))

    result <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
    result_sv <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))

    for (i in 1:nrow(data1)) {
        pos_ind <- which(data1_pos_uniq_sort[, 1] == data1[i, 1] & data1_pos_uniq_sort[, 2] == data1[i, 2] & data1_pos_uniq_sort[, 3] == data1[i, 3])
        cell_ind <- which(data1_cell_uniq_sort[, 1] == data1[i, 5])
        result[cell_ind, pos_ind] <- data1[i, 13]
        result_sv[cell_ind, pos_ind] <- data1[i, 15]
    }
    rownames(result) <- data1_cell_uniq_sort
    colnames(result) <- data1_pos_uniq_sort$posind
    result[result == Inf] <- max(c(1, result[result != Inf]))
    rownames(result_sv) <- data1_cell_uniq_sort
    colnames(result_sv) <- data1_pos_uniq_sort$posind

    breaksList <- seq(4, 30, by = 0.1)
    breaksList <- append(breaksList, max(result))
    breaksList <- append(breaksList, -1, 0)
    breaksList <- unique(sort(breaksList))
    mycol <- colorpanel(n = length(breaksList) - 1, low = "lightgoldenrodyellow", mid = "darkorange", high = "firebrick2")

    res <- pheatmap(result, border_color = NA, show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = T, clustering_method = "ward.D", scale = "none", col = mycol, cex = 0.5, main = inputfile, breaks = breaksList)

    ## Assign sv calls to the genome-wide bins (200kb)
    par("mar")
    par(mar = c(0.5, 0.5, 0.5, 0.5))

    genome_bins_sort <- genome_bins[genome_bins[, 1] == chrom[1], ]
    for (i in 2:length(chrom)) {
        genome_bins_sort <- rbind(genome_bins_sort, genome_bins[genome_bins[, 1] == chrom[i], ])
    }
    print(genome_bins_sort)
    print(result_sv)
    print(data1_pos_uniq_sort)
    result_sv_bins <- matrix(0, nrow(genome_bins_sort), nrow(result_sv))
    for (i in 1:ncol(result_sv)) {
        num_bins1 <- which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 2] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 2])
        print(num_bins1)
        num_bins2 <- which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 3] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 3])
        print(num_bins2)
        if (sum(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 3] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 3]) == 0) {
            num_bins2 <- max(which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1])))
        }

        num_bins <- c(num_bins1:num_bins2)
        for (j in 1:length(num_bins)) {
            result_sv_bins[num_bins[j], ] <- t(as.matrix(result_sv[, i]))
        }
        cat(paste0(i, " "))
    }

    widths <- matrix(0, length(chrom), 1)
    for (i in 1:length(chrom)) {
        widths[i] <- sum(genome_bins_sort[, 1] == chrom[i]) / 100
    }

    colnames(genome_bins_sort) <- colnames(data1_pos_uniq_sort)
    gos2 <- genome_bins_sort
    mat2 <- result_sv_bins[, res$tree_row$order]
    mat2 <- mat2[, ncol(mat2):1]
    dnull <- matrix(0, nrow(mat2), ncol(mat2))

    l <- layout(matrix(seq(1, 23), 1, 23, byrow = TRUE), widths = widths)
    for (i in 1:length(chrom)) {
        d <- mat2[gos2$chrom == chrom[i], ]
        if (sum(d) == 0) {
            image(seq_len(nrow(dnull)), seq_len(ncol(dnull)), dnull, zlim = c(0, 12), col = "white", xlab = "", ylab = "", axes = FALSE, main = chrom[i], cex.main = 0.8)
            box()
        }
        if (sum(d) > 1) {
            image(seq_len(nrow(d)), seq_len(ncol(d)), d, zlim = c(0, 12), col = colors, xlab = "", ylab = "", axes = FALSE, main = chrom[i], cex.main = 0.8)
            box()
        }
    }

    dev.off()
}