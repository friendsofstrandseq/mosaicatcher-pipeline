


# plot.clustering <- function(inputfile, bin.bed.filename, position.outputfile, chromosome.outputfile,) {
# plot.clustering <- function(inputfile, bin.bed.filename, position.outputfile, chromosome.outputfile, chromosomes) {
plot.clustering <- function(inputfile, bin.bed.filename, position.outputfile, chromosomes) {
 

    library(pheatmap)
    library(gplots)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(reshape2)


    # Colors, SV types & Chroms instanciation
    ash12rainbow <-
      c(
        "#77AADD",
        "#4477AA",
        "#114477",
        "#CC99BB",
        "#AA4488",
        "#771155",
        "#DDDD77",
        "#AAAA44",
        "#777711",
        "#DDAA77",
        "#AA7744",
        "#774411"
      )
    sv_call_name <-
      c(
        "del_h1",
        "del_h2",
        "del_hom",
        "dup_h1",
        "dup_h2",
        "dup_hom",
        "inv_h1",
        "inv_h2",
        "inv_hom",
        "idup_h1",
        "idup_h2",
        "complex"
      )
    chrom <-
      c(
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX"
      )
    
    chrom <- chromosomes
    print(chrom)
    # stop()
    
    # Genome bins loading
    ## Correspond to 200kb size bins for each of the chroms
    genome_bins <-
      read.table(
        bin.bed.filename,
        sep = "\t",
        header = F,
        comment.char = ""
      )
    
    
    
    # READ SV calls data
    data1 <-
      read.table(inputfile,
                 sep = "\t",
                 header = T,
                 comment.char = "")
    
    
    data1$color <- 0
    
    
    
    # GOBACK
    list_directory <- dir("./", full.names = TRUE)
    
    # Open PDF for output file
    pdf(position.outputfile, width = 11, height = 10)
    
    # FOR loop to replace AF by the corresponding index in sv_call_name list
    for (j in 1:length(sv_call_name)) {
      # Select all elements corresponding to SV index
      tmp <- which(data1[, 9] == sv_call_name[j])
      # AF column value become index list
      data1[tmp, 15] <- j
    }
    
    # Subset only SV segments positions (a SV can overlap multiple segments)
    data1_pos <- data1[, 1:3]
    # Unique
    data1_pos_uniq <- unique(data1_pos)
    # Sort only chr1 (numeric sorting)
    data1_pos_uniq_sort <-
      data1_pos_uniq[data1_pos_uniq$chrom == "chr1", ]
    
    # Do the same for chrom 2 to X & concatenate at each iteration
    for (i in 2:length(chrom)) {
      data1_pos_uniq_sort <-
        rbind(data1_pos_uniq_sort, data1_pos_uniq[data1_pos_uniq$chrom == chrom[i],])
    }
    # New column to get position in the list
    data1_pos_uniq_sort$posind <- c(1:nrow(data1_pos_uniq_sort))
    
    # Get all cell names
    data1_cell <- data1[, 5]
    # Unique
    data1_cell_uniq <- unique(data1_cell)
    # Sort
    data1_cell_uniq_sort <- as.matrix(sort(data1_cell_uniq))
    
    # Matrix, row=nb of cells, col=nb of SV (unique positions)
    result <-
      matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
    
    # Matrix, row=nb of cells, col=nb of SV (unique positions)
    result_sv <-
      matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
    
    print(result)
    print(result_sv)
    print(data1)
    
    # FOR loop: for each SV
    for (i in 1:nrow(data1)) {
      # Get pos index in data1_pos_uniq_sort
      pos_ind <-
        which(
          data1_pos_uniq_sort[, 1] == data1[i, 1] &
            data1_pos_uniq_sort[, 2] == data1[i, 2] &
            data1_pos_uniq_sort[, 3] == data1[i, 3]
        )
      # Get cell index in data1_pos_uniq_sort
      cell_ind <- which(data1_cell_uniq_sort[, 1] == data1[i, 5])
      

      # Fill the matrix at corresponding indexes with llr_to_ref & af (SV type) columns
      result[cell_ind, pos_ind] <- data1[i, 13]
      result_sv[cell_ind, pos_ind] <- data1[i, 15]
    }
    
    # Replace row & col names
    rownames(result) <- data1_cell_uniq_sort
    colnames(result) <- data1_pos_uniq_sort$posind
    rownames(result_sv) <- data1_cell_uniq_sort
    colnames(result_sv) <- data1_pos_uniq_sort$posind
    
    # Replace Inf values
    result[result == Inf] <- max(c(1, result[result != Inf]))
    
    
    
    
    ########
    
    
    # Create seq
    breaksList <- seq(4, 30, by = 0.1)
    # Add max of result matrix
    breaksList <- append(breaksList, max(result))
    # Add -1 at index 0
    breaksList <- append(breaksList,-1, 0)
    # Unique & Sort
    breaksList <- unique(sort(breaksList))
    # Create color range
    mycol <-
      gplots::colorpanel(
        n = length(breaksList) - 1,
        low = "lightgoldenrodyellow",
        mid = "darkorange",
        high = "firebrick2"
      )
    
    
    # Sed seed for R Color Brewer
    set.seed(90)
    
    
    # Generate color list based on chroms length
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    colors_chroms = sample(col_vector, length(chrom))
    # Instanciate data$color column
    data1_pos_uniq_sort$color = "NA"
    # Iterate over chrom to attribute color
    for (i in 1:length(chrom)) {
      chrom_index_list <- which(data1_pos_uniq_sort$chrom == chrom[i])
      data1_pos_uniq_sort[chrom_index_list, "color"] = colors_chroms[i]
    }
    
    # Sort data based on define factor
    chrOrder <-
      c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
    data1_pos_uniq_sort$chrom <-
      factor(data1_pos_uniq_sort$chr, levels = chrOrder)
    data1_pos_uniq_sort[order(data1_pos_uniq_sort$chrom),]
    
    # Retrieve information & metadata for ComplexHeatmap
    Chroms <- sample(col_vector, length(chrom))
    names(Chroms) <- chrom
    anno_colors <- list(Chroms = Chroms)
    col_annotation <- as.data.frame(data1_pos_uniq_sort$chrom)
    rownames(col_annotation) <- colnames(result)
    colnames(col_annotation) <- "Chroms"
    
    # Plot Clustered Heatmap
    res <-
      pheatmap::pheatmap(
        # name = "HELLO",
        result,
        border_color = NA,
        show_rownames = T,
        show_colnames = F,
        cluster_cols = F,
        cluster_rows = T,
        clustering_method = "ward.D",
        scale = "none",
        col = mycol,
        main = tools::file_path_sans_ext(basename(inputfile)),
        annotation_col = col_annotation,
        annotation_colors = anno_colors,
        breaks = breaksList,
        annotation_legend = TRUE
      )
    # print(res)
    # dev.off()
    
    ######################
    
    # Add chrom name + start_end
    chr_name <- matrix("chr", ncol(result), 1)
    for (i in 1:ncol(result)) {
      chr_name[i, 1] <-
        paste0(data1_pos_uniq_sort[i, 1],
               "_",
               data1_pos_uniq_sort[i, 2],
               "_",
               data1_pos_uniq_sort[i, 3])
    }

    # Add as col names to result_sv    
    colnames(result_sv) <- chr_name
    
    # SV list
    sv_list <-
      c(
        "none",
        "del_h1",
        "del_h2",
        "del_hom",
        "dup_h1",
        "dup_h2",
        "dup_hom",
        "inv_h1",
        "inv_h2",
        "inv_hom",
        "idup_h1",
        "idup_h2",
        "complex"
      )

    # SV type colors
    colors <-
      structure(
        c(
          "white",
          "#77AADD",
          "#4477AA",
          "#114477",
          "#CC99BB",
          "#AA4488",
          "#771155",
          "#DDDD77",
          "#AAAA44",
          "#777711",
          "#DDAA77",
          "#AA7744",
          "#774411"
        ),
        names = c(0:12)
      )
    
 
    # Create data frame based on previous pheatmap
    mat <- as.data.frame(result_sv[res$tree_row$order, ])
    # library(reshape2)
    melt_mat <- reshape2::melt(data1$af)
    print(sv_call_name)
    print(data1$af)
    print(data1$sv_call_name)
    print(table(melt_mat$value))
    sv_list_sub = sv_list[as.numeric(names(table(melt_mat$value)))]
    print(sv_list_sub)
    new_colors = colors[as.numeric(names(table(melt_mat$value)))]
    
    # 
    # 
    # print(head(mat))
    # print(colors)
    # print(sv_list_sub)
    # print(length(sv_list_sub))
    # print(new_colors)
    # print(length(new_colors))
    # # print(mat)
    # print(1:length(sv_list_sub))

    print(colors)
    print(typeof(colors))
    print(new_colors)
    print(sv_list_sub)
    print(length(sv_list_sub))
    # print(0:length(sv_list_sub)-1)
    len_sv_list_sub = length(sv_list_sub) - 1
    print(len_sv_list_sub)
    print(0:len_sv_list_sub)
    print(length(0:len_sv_list_sub))
    # Heatmap metadata
    ha_column <-
      ComplexHeatmap::HeatmapAnnotation(df = col_annotation,
                                        col = anno_colors)
    print(mat)
    print(length(sv_list))
    print(0:(length(sv_list)-1))

    # Draw Heatmap based on SV type
    ht1 <-
      ComplexHeatmap::Heatmap(
        data.matrix(mat),
        name = "SV type",
        col = colors,
        heatmap_legend_param = list(at = 0:(length(sv_list)-1), labels = sv_list),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = tools::file_path_sans_ext(basename(inputfile)),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 1),
        column_title_gp = gpar(fontsize = 7, fontface = "bold"),
        top_annotation = ha_column
      )
    # print(ht1)
    draw(ht1, show_annotation_legend = TRUE)
     

    # FIXME: second part of plots are not functional

    # dev.off()

    # pdf(chromosome.outputfile, width = 11, height = 5)

    # data1 <- read.table(inputfile, sep = "\t", header = T, comment.char = "")

    # data1$color <- 0
    # ash12rainbow <- c("#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411")
    # sv_call_name <- c("del_h1", "del_h2", "del_hom", "dup_h1", "dup_h2", "dup_hom", "inv_h1", "inv_h2", "inv_hom", "idup_h1", "idup_h2", "complex")
    # for (j in 1:length(sv_call_name)) {
    #     tmp <- which(data1[, 9] == sv_call_name[j])
    #     data1[tmp, 15] <- j
    # }
    # colors <- structure(c("white", "#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411"), names = c(0:12))


    # data1_pos <- data1[, 1:3]
    # data1_pos_uniq <- unique(data1_pos)
    # data1_pos_uniq_sort <- data1_pos_uniq[data1_pos_uniq$chrom == "chr1", ]
    # # chrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

    # for (i in 2:length(chrom)) {
    #     data1_pos_uniq_sort <- rbind(data1_pos_uniq_sort, data1_pos_uniq[data1_pos_uniq$chrom == chrom[i], ])
    # }

    # data1_pos_uniq_sort$posind <- c(1:nrow(data1_pos_uniq_sort))

    # data1_cell <- data1[, 5]
    # data1_cell_uniq <- unique(data1_cell)
    # data1_cell_uniq_sort <- as.matrix(sort(data1_cell_uniq))

    # result <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))
    # result_sv <- matrix(0, nrow(data1_cell_uniq_sort), nrow(data1_pos_uniq))

    # for (i in 1:nrow(data1)) {
    #     pos_ind <- which(data1_pos_uniq_sort[, 1] == data1[i, 1] & data1_pos_uniq_sort[, 2] == data1[i, 2] & data1_pos_uniq_sort[, 3] == data1[i, 3])
    #     cell_ind <- which(data1_cell_uniq_sort[, 1] == data1[i, 5])
    #     result[cell_ind, pos_ind] <- data1[i, 13]
    #     result_sv[cell_ind, pos_ind] <- data1[i, 15]
    # }
    # rownames(result) <- data1_cell_uniq_sort
    # colnames(result) <- data1_pos_uniq_sort$posind
    # result[result == Inf] <- max(c(1, result[result != Inf]))
    # rownames(result_sv) <- data1_cell_uniq_sort
    # colnames(result_sv) <- data1_pos_uniq_sort$posind

    # breaksList <- seq(4, 30, by = 0.1)
    # breaksList <- append(breaksList, max(result))
    # breaksList <- append(breaksList, -1, 0)
    # breaksList <- unique(sort(breaksList))
    # mycol <- colorpanel(n = length(breaksList) - 1, low = "lightgoldenrodyellow", mid = "darkorange", high = "firebrick2")

    # res <- pheatmap(result, border_color = NA, show_rownames = T, show_colnames = F, cluster_cols = F, cluster_rows = T, clustering_method = "ward.D", scale = "none", col = mycol, cex = 0.5, main = inputfile, breaks = breaksList)

    # ## Assign sv calls to the genome-wide bins (200kb)
    # par("mar")
    # par(mar = c(0.5, 0.5, 0.5, 0.5))

    # genome_bins_sort <- genome_bins[genome_bins[, 1] == chrom[1], ]
    # for (i in 2:length(chrom)) {
    #     genome_bins_sort <- rbind(genome_bins_sort, genome_bins[genome_bins[, 1] == chrom[i], ])
    # }
    # print(genome_bins_sort)
    # print(result_sv)
    # print(data1_pos_uniq_sort)
    # result_sv_bins <- matrix(0, nrow(genome_bins_sort), nrow(result_sv))
    # for (i in 1:ncol(result_sv)) {
    #     num_bins1 <- which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 2] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 2])
    #     print(num_bins1)
    #     num_bins2 <- which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 3] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 3])
    #     print(num_bins2)
    #     if (sum(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1]) & genome_bins_sort[, 2] <= data1_pos_uniq_sort[i, 3] & genome_bins_sort[, 3] > data1_pos_uniq_sort[i, 3]) == 0) {
    #         num_bins2 <- max(which(genome_bins_sort[, 1] == as.character(data1_pos_uniq_sort[i, 1])))
    #     }

    #     num_bins <- c(num_bins1:num_bins2)
    #     for (j in 1:length(num_bins)) {
    #         result_sv_bins[num_bins[j], ] <- t(as.matrix(result_sv[, i]))
    #     }
    #     cat(paste0(i, " "))
    # }

    # widths <- matrix(0, length(chrom), 1)
    # for (i in 1:length(chrom)) {
    #     widths[i] <- sum(genome_bins_sort[, 1] == chrom[i]) / 100
    # }

    # colnames(genome_bins_sort) <- colnames(data1_pos_uniq_sort)
    # gos2 <- genome_bins_sort
    # mat2 <- result_sv_bins[, res$tree_row$order]
    # mat2 <- mat2[, ncol(mat2):1]
    # dnull <- matrix(0, nrow(mat2), ncol(mat2))

    # l <- layout(matrix(seq(1, 23), 1, 23, byrow = TRUE), widths = widths)
    # for (i in 1:length(chrom)) {
    #     d <- mat2[gos2$chrom == chrom[i], ]
    #     if (sum(d) == 0) {
    #         image(seq_len(nrow(dnull)), seq_len(ncol(dnull)), dnull, zlim = c(0, 12), col = "white", xlab = "", ylab = "", axes = FALSE, main = chrom[i], cex.main = 0.8)
    #         box()
    #     }
    #     if (sum(d) > 1) {
    #         image(seq_len(nrow(d)), seq_len(ncol(d)), d, zlim = c(0, 12), col = colors, xlab = "", ylab = "", axes = FALSE, main = chrom[i], cex.main = 0.8)
    #         box()
    #     }
    # }

    dev.off()
}