suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

plotNBdist <- function(probs, alpha = 0.05) {
  
  assert_that(is.data.table(probs),
              "sample" %in% colnames(probs),
              "cell"   %in% colnames(probs),
              "chrom"  %in% colnames(probs),
              "start"  %in% colnames(probs),
              "end"    %in% colnames(probs),
              "class"  %in% colnames(probs),
              "nb_p"   %in% colnames(probs),
              "expected" %in% colnames(probs),
              "W"      %in% colnames(probs),
              "C"      %in% colnames(probs)  )
  
  assert_that(length(unique(probs$cell)) <= 6,
              nrow(unique(probs[, .(chrom, start, end)])) <= 6)
  
  # make a copy in order to not overwrite the original table
  x <- copy(probs)
  
  # reduce table to 1 entry per locus and cell
  x <- unique(x[, .(sample, cell, chrom, start, end, class, nb_p, expected, W, C)])
  x[, assert_that(.N == 1), by = .(chrom, start, end, sample, cell)] %>% invisible

  # add nb_r column
  x[, nb_r := expected * nb_p / (1-nb_p)]
  
  # expand each cell and segment by a range to plot the NB distribution
  x[, `:=`(nb_range_start = 0,
           nb_range_step  = round(log(expected)),
           nb_range_end   = round(2*expected)) ]
  x <- x[, 
         .(nb_x = seq(nb_range_start, nb_range_end, nb_range_step),
           nb_cn0 = dnbinom(seq(nb_range_start, nb_range_end, nb_range_step), size = alpha * nb_r, prob = nb_p), 
           nb_cn1 = dnbinom(seq(nb_range_start, nb_range_end, nb_range_step), size = 0.5 * nb_r, prob = nb_p),
           nb_cn2 = dnbinom(seq(nb_range_start, nb_range_end, nb_range_step), size = nb_r, prob = nb_p)),
         by = .(chrom, start, end, sample, cell, expected, W, C)]
  
  plt <- ggplot(x) + 
    geom_line(aes(nb_x, nb_cn0), color = "red") + 
    geom_line(aes(nb_x, nb_cn1), color = "green") + 
    geom_line(aes(nb_x, nb_cn2), color = "blue") + 
    geom_vline(aes(xintercept = W), linetype = "dashed", color = "sandybrown") +
    geom_vline(aes(xintercept = C), linetype = "dashed", color = "paleturquoise4") +
    facet_grid(chrom + start + end ~ sample + cell, scales = "free")
    
  return(plt)
}




# test case:
# x <- probs[cell %in% c("cell_0","cell_1","cell_2","cell_3") & chrom=="chr10" & start < 30e6]
# plt <- plotNBdist(x)
# cairo_pdf("test.pdf", width = 20, height = 16)
# print(plt)
# dev.off()
