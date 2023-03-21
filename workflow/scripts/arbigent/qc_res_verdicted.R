# This file QCs the res_verdicted, e.g. to see
# problems with normalization etc.

library(stringr)
library(ggplot2)
library(optparse)
library(reshape2)
library(tibble)
library(matrixStats)
library(grid)
library(dplyr)

make_barplot_invwise <- function(rv_f, samples_f) {
  number <- 1
  ### Define events again ###
  ref_events <- c("0/0", "0|0")
  simple_events <- c("0|1", "0/1", "1|0", "1/0", "1|1", "1/1")
  simple_events_lowconf <- c(
    "0/0_lowconf", "0|0_lowconf", "0|1_lowconf",
    "0/1_lowconf", "1|0_lowconf", "1|0_lowconf", "1/0_lowconf", "1|1_lowconf", "1/1_lowconf"
  )
  verdicts <- c("simple", "simple_lowconf", "ref", "noreads", "complex")

  ### Replace GTs with 'event' - simple, ref, complex, noreads... ###
  res_verdict <- rv_f %>%
    mutate_all(funs(ifelse(. %in% simple_events, "simple", .))) %>%
    mutate_all(funs(ifelse(. %in% simple_events_lowconf, "simple_lowconf", .))) %>%
    mutate_all(funs(ifelse(. %in% ref_events, "ref", .))) %>%
    mutate_all(funs(ifelse(. %in% c("simple", "simple_lowconf", "ref", "noreads"), ., "complex")))

  ### We accidentally replaced chr, start, end too, so we want to recover them.
  res_verdict$chrom <- rv_f$chrom
  res_verdict$start <- rv_f$start
  res_verdict$end <- rv_f$end
  res_verdict$verdict <- rv_f$verdict

  # Now count which verdict is how frequent per inversion. We need this only for one
  # line (marked with !!!##!!!)
  for (verdict in verdicts) {
    res_verdict[[verdict]] <- rowCounts(as.matrix(res_verdict), value = verdict)
  }

  # Enumerate inversions
  res_verdict$invno <- row.names(res_verdict)

  res_verdict_verdictsorted <- res_verdict[order(res_verdict$verdict, res_verdict$simple), "invno"]
  # This is that line: !!!##!!!
  # Get the sorted verdicts (this is a vector of length n_inversions)
  verd <- res_verdict[res_verdict_verdictsorted, "verdict"]
  # Get the position at which we jump to the next verdict
  xpos <- c(0, cumsum(rle(verd)$lengths))
  # Get all verdictnames
  xnames <- rle(verd)$values

  # We go towards plotting. #
  baseheight <- length(samples_f)
  n_classes <- length(xnames)


  res_verdict_molten <- melt(res_verdict[, c("invno", "verdict", samples), drop = F], id.vars = c("invno", "verdict"))
  print("got so far")
  ### PLOT ###
  p2 <- ggplot() +
    geom_bar(data = res_verdict_molten, aes(x = as.character(invno), fill = value)) +
    scale_fill_manual(values = c("blue", "white", "darkgrey", "darkgreen", "green")) +
    xlim(res_verdict_verdictsorted)

  p3 <- p2 + geom_segment(
    aes(
      y = baseheight,
      yend = baseheight,
      x = xpos[1:length(xpos) - 1],
      xend = xpos[2:length(xpos)],
      color = xnames[1:length(xnames)]
    ),
    arrow = arrow(ends = "both"),
  ) +
    geom_text(
      aes(
        x = ((xpos[2:length(xpos)] - xpos[1:length(xpos) - 1]) / 2) + xpos[1:length(xpos) - 1],
        y = seq(from = baseheight + 2, to = baseheight + 2 + (n_classes * 1.3), length.out = n_classes),
        label = xnames[1:length(xnames)],
        color = xnames[1:length(xnames)]
      ),
      check_overlap = F,
      size = 3,
      fontface = "bold",
      angle = 0
    )

  return(p3)
}


make_barplot_groupwise <- function(idf_f) {
  xorder <- as.character((idf_f[order(idf_f$simpleGT, decreasing = T), ]$variable))
  idf_f <- within(idf_f, variable <- factor(variable, levels = xorder))

  idf_2 <- idf_f %>%
    group_by(simpleGT) %>%
    mutate(sum_events = sum(value)) %>%
    slice(1)
  p2 <- ggplot(idf_2) +
    geom_bar(aes(x = simpleGT, y = sum_events, fill = simpleGT), stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(rev(xorder))

  return(p2)
}

make_barplot_all <- function(idf_f) {
  xorder <- as.character((idf_f[order(idf_f$simpleGT, decreasing = T), ]$variable))

  idf_f <- within(idf_f, variable <- factor(variable, levels = xorder))

  p1 <- ggplot(idf_f) +
    geom_bar(aes(x = variable, y = value, fill = simpleGT), stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(rev(xorder))

  return(p1)
}


make_inventory <- function(rv_f, samples_f) {
  # invent = inventory. We flatten the matrix and count elements
  invent <- table(unlist(rv_f[, samples_f, drop = F]))
  idf <- melt(as.data.frame.matrix(rbind((invent))))

  # What to color how?
  idf$simpleGT <- "complex"
  ref_events <- c("0/0", "0|0")
  simple_events <- c("0|1", "0/1", "1|0", "1/0", "1|1", "1/1")
  simple_events_lowconf <- c(
    "0/0_lowconf", "0|0_lowconf", "0|1_lowconf",
    "0/1_lowconf", "1|0_lowconf", "1|0_lowconf", "1/0_lowconf", "1|1_lowconf", "1/1_lowconf"
  )

  idf[idf$variable %in% ref_events, "simpleGT"] <- "a) ref"
  idf[idf$variable %in% simple_events, "simpleGT"] <- "b) inv_simple"
  idf[idf$variable %in% simple_events_lowconf, "simpleGT"] <- "c) inv_simple_lowconf"
  idf[idf$variable == "noreads", "simpleGT"] <- "z) noreads"

  return(idf)
}




# INPUT INSTRUCTIONS
option_list <- list(
  make_option(c("-f", "--file"),
    type = "character", default = NULL,
    help = "res_verdicted", metavar = "character"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./outputcorr/",
    help = "Outputdir for phased all.txt and other qcs", metavar = "character"
  )
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
res_verdicted_link <- opt$file
outdir <- opt$outdir

# res_verdicted_link = "~/s/g/korbel/hoeps/projects/huminvs/mosai_results/results_freeze4manual/regenotyper_allsamples_bulk/arbigent_results/res_verdicted.vcf"
# res_verdicted_link = "/home/hoeps/Desktop/hufsah_freeze4manual/Arbigent_gts.vcf"
# outdir = '~/Desktop/hufsah_freeze4manual/'

############# RUN CODE #################

# Load res_verdicted
rv <- read.table(res_verdicted_link, header = 1, sep = "\t", stringsAsFactors = F)

# Get the samplenames
# samples <- colnames(rv)[grep("^[HMNG].*", colnames(rv))]
samples <- tail(colnames(rv), 1)

print(rv)
print(samples)
# First, inventory.
idf <- make_inventory(rv, samples)

print(idf)
# Give numbers
gtclass_merge <- (idf %>% group_by(simpleGT) %>% mutate(sum_events = sum(value)) %>% slice(1))[, c("simpleGT", "sum_events")]
verdicts_invwise <- as.data.frame.matrix(rbind(table(rv$verdict)))

# Make plots
p1 <- make_barplot_all(idf)
p2 <- make_barplot_groupwise(idf)
p3 <- make_barplot_invwise(rv, samples)

### SAVE ###

# a) tables
write.table(gtclass_merge, file = paste0(outdir, "/simple-complex-numbers.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(t(verdicts_invwise), file = paste0(outdir, "/verdicts-numbers.txt"), sep = "\t", row.names = T, col.names = F, quote = F)

# b) plots
ggsave(file = paste0(outdir, "/all_gts_overview.pdf"), plot = p1, width = 30, height = 10, units = "cm", device = "pdf")
ggsave(file = paste0(outdir, "/all_verdicts_overview.pdf"), plot = p2, width = 15, height = 10, units = "cm", device = "pdf")
ggsave(file = paste0(outdir, "/lineplot_gts.pdf"), plot = p3, width = 30, height = 15, units = "cm", device = "pdf")
