#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
# options(error = function() traceback(3))


suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(scales) %>% invisible()
library(assertthat) %>% invisible()
library(stringr) %>% invisible()
library(RColorBrewer) %>% invisible()

################################################################################
# Settings                                                                     #
################################################################################

zcat_command <- "zcat"
# FIXME : tmp solution
chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

format_Mb <- function(x) {
  paste(comma(x / 1e6), "Mb")
}

### Colors for background
manual_colors <- c(
  # duplications
  simul_hom_dup = "firebrick4",
  dup_hom = muted("firebrick4", 70, 50),
  simul_het_dup = "firebrick2",
  dup_h1 = muted("firebrick2", 90, 30),
  dup_h2 = muted("firebrick2", 80, 20),
  # deletions
  simul_hom_del = "dodgerblue4",
  del_hom = muted("dodgerblue4", 50, 60),
  simul_het_del = "dodgerblue2",
  del_h1 = muted("dodgerblue2", 80, 50),
  del_h2 = muted("deepskyblue2", 80, 50),
  # inversions
  simul_hom_inv = "chartreuse4",
  inv_hom = muted("chartreuse4", 80, 50),
  simul_het_inv = "chartreuse2",
  inv_h1 = muted("chartreuse2", 100, 60),
  inv_h2 = muted("darkolivegreen3", 100, 60),
  # other SVs
  simul_false_del = "darkgrey",
  simul_inv_dup = "darkgoldenrod2",
  idup_h1 = muted("darkgoldenrod2", 80, 70),
  idup_h2 = muted("gold", 80, 70),
  complex = "darkorchid1",
  # background
  bg1 = "#ffffff",
  bg2 = "#aaafaa",
  # Strand states
  `State: WW` = "sandybrown",
  `State: CC` = "paleturquoise4",
  `State: WC` = "khaki",
  `State: CW` = "yellow2"
)




################################################################################
# Usage                                                                        #
################################################################################

print_usage_and_stop <- function(msg = NULL) {
  if (!is.null(msg)) {
    message(msg)
  }
  message("Plot Strand-seq counts of all cells for a single chromosome.                    ")
  message("                                                                                ")
  message("Usage:                                                                          ")
  message("    Rscript chrom.R [OPTIONS] <count-file> <chrom> <out.pdf>                    ")
  message("                                                                                ")
  message("OPTIONS (no spaces around `=`):                                                 ")
  message("    per-page=<int>            Number of cells to be printed per page            ")
  message("    segments=<file>           Show the segmentation in the plots                ")
  message("    singlecellsegments=<file> Show per-cell  segmentation in the plots          ")
  message("    calls=<file>              Highlight SV calls provided in a table            ")
  message("    truth=<file>              Mark the `true`` SVs provided in a table          ")
  message("    strand=<file>             Mark the strand states which calls are based on   ")
  message("    complex=<file>            Mark complex regions given in file                ")
  message("    groups=<file>             Table with SV call grouping                       ")
  message("    no-none                   Do not hightlight black-listed (i.e. None) bins   ")
  message("                                                                                ")
  message("Generates one plot per chromosome listing all cells below another, separated    ")
  message("into pages. If an SV probability file is provided (2), segments are colored     ")
  message("according to the predicted SV classes. Note that only certain classes are       ")
  message("accepted and the script will yield an error if others are provided.             ")
  message("Similarly, a segmentation file can be specified, yet it must contain exactly    ")
  message("one segmentation (MosaiCatcher reports segmentations for various total numbers  ")
  message("of breakpoints).")
  options(show.error.messages = F)
  stop()
}



################################################################################
# Command Line Arguments                                                       #
################################################################################

args <- commandArgs(trailingOnly = T)

if (length(args) < 3) print_usage_and_stop("[Error] Too few arguments!")

f_counts <- args[length(args) - 2]
CHROM <- args[length(args) - 1]
f_out <- args[length(args)]

f_segments <- NULL
f_calls <- NULL
f_truth <- NULL
f_strand <- NULL
f_complex <- NULL
cells_per_page <- 8
show_none <- T

if (length(args) > 3) {
  if (!all(grepl("^(strand|calls|segments|per-page|truth|no-none|complex|singlecellsegments|groups)=?", args[1:(length(args) - 3)]))) {
    print_usage_and_stop("[Error]: Options must be one of `calls`, `segments`, `per-page`, or `truth`")
  }
  for (op in args[1:(length(args) - 3)]) {
    if (grepl("^segments=", op)) f_segments <- str_sub(op, 10)
    if (grepl("^calls=", op)) f_calls <- str_sub(op, 7)
    if (grepl("^truth=", op)) f_truth <- str_sub(op, 7)
    if (grepl("^per-page=", op)) {
      pp <- as.integer(str_sub(op, 10))
      if (pp > 0 && pp < 50) {
        cells_per_page <- pp
      }
    }
    if (grepl("^strand=", op)) f_strand <- str_sub(op, 8)
    if (grepl("^complex=", op)) f_complex <- str_sub(op, 9)
    if (grepl("^groups=", op)) f_groups <- str_sub(op, 8)
    if (grepl("^singlecellsegments=", op)) f_scsegments <- str_sub(op, 20)
    if (grepl("^no-none$", op)) show_none <- F
  }
}


################################################################################
# Read & check input data                                                      #
################################################################################

### Check counts table
message(" * Reading count data ", f_counts, "...")
if (grepl("\\.gz$", f_counts)) {
  counts <- fread(paste(zcat_command, f_counts))
} else {
  counts <- fread(f_counts)
}

# FIXME : tmp
# print(counts)
counts <- counts[counts$chrom %in% chroms, ]
# print(counts)

assert_that(
  "chrom" %in% colnames(counts),
  "start" %in% colnames(counts),
  "end" %in% colnames(counts),
  "class" %in% colnames(counts),
  "sample" %in% colnames(counts),
  "cell" %in% colnames(counts),
  "w" %in% colnames(counts),
  "c" %in% colnames(counts)
) %>% invisible()
counts[, sample_cell := paste(sample, "-", cell)]
setkey(counts, chrom, sample_cell)
bins <- unique(counts[, .(chrom, start, end)])

### Check CHROM:
assert_that(CHROM %in% unique(counts$chrom)) %>% invisible()
counts <- counts[chrom == CHROM]




### Check SV call file
if (!is.null(f_calls)) {
  message(" * Reading SV calls from ", f_calls, "...")
  svs <- fread(f_calls)
  assert_that(
    "chrom" %in% colnames(svs),
    "start" %in% colnames(svs),
    "end" %in% colnames(svs),
    "sample" %in% colnames(svs),
    "cell" %in% colnames(svs),
    ("SV_class" %in% colnames(svs) | "sv_call_name" %in% colnames(svs))
  ) %>% invisible()
  if (!("SV_class" %in% colnames(svs))) {
    svs[, SV_class := sv_call_name]
  }
  assert_that(all(svs$SV_class %in% names(manual_colors))) %>% invisible()
  svs[, sample_cell := paste(sample, "-", cell)]

  set_diff <- setdiff(unique(svs$sample_cell), unique(counts$sample_cell))
  if (length(set_diff) > 0) message("[Warning] SV calls and Counts differ in cells: ", set_diff)

  svs <- svs[chrom == CHROM]
}

### Check segment table
if (!is.null(f_segments)) {
  message(" * Reading segmentation file from ", f_segments, "...")
  seg <- fread(f_segments)


  # FIXME : tmp
  seg <- seg[seg$chrom %in% chroms, ]
  seg <- seg[seg$bps > 0, ] # SOLVE T2T ISSUE

  # print(seg)

  assert_that(
    "chrom" %in% colnames(seg),
    "bps" %in% colnames(seg)
  ) %>% invisible()
  if ("k" %in% colnames(seg)) {
    seg[, assert_that(length(unique(k)) == 1), by = .(chrom)] %>% invisible()
  }

  # print(seg)

  # print(bins)
  # print(bins[, .N, by = chrom])
  # # print(bins %>% count(chrom))
  # print(bins[, .N, by = chrom][, .(chrom, N = c(0, cumsum(N)))])
  # print(bins[, .N, by = chrom][, .(chrom, N = c(0, cumsum(N))[1:(.N - 1)])])



  seg <- merge(seg, bins[, .N, by = chrom][, .(chrom, N = c(0, cumsum(N))[1:(.N - 1)])], by = "chrom")

  # print(c(1, bps[1:(.N - 1)] + 1))
  # print(bps)
  # print(seg)

  # stop()


  seg[, `:=`(from = c(1, bps[1:(.N - 1)] + 1), to = bps), by = chrom]

  # print(seg)


  seg[, `:=`(
    start = bins[from + N]$start,
    end = bins[to + N]$end
  )]

  # print("TEST")


  seg[, SV_class := rep(c("bg1", "bg2"), .N)[1:.N], by = chrom]



  seg <- seg[chrom == CHROM]
}



### Check simulated variants
if (!is.null(f_truth)) {
  message(" * Reading simulated variants from ", f_truth, "...")
  simul <- fread(f_truth)
  assert_that(
    "chrom" %in% colnames(simul),
    "start" %in% colnames(simul),
    "end" %in% colnames(simul),
    "sample" %in% colnames(simul),
    "cell" %in% colnames(simul),
    "SV_type" %in% colnames(simul)
  ) %>% invisible()
  simul[, `:=`(SV_class = paste0("simul_", SV_type), SV_type = NULL, sample_cell = paste(sample, "-", cell))]
  simul[, sample_cell := paste(sample, "-", cell)]

  set_diff <- setdiff(unique(simul$sample_cell), unique(counts$sample_cell))
  if (length(set_diff) > 0) message("[Warning] True SVs and Counts differ in cells: ", set_diff)

  simul <- simul[chrom == CHROM]
}

### Check strand states file
if (!is.null(f_strand)) {
  message(" * Reading strand state file from ", f_strand, "...")
  strand <- fread(f_strand)
  assert_that(
    "sample" %in% colnames(strand),
    "cell" %in% colnames(strand),
    "chrom" %in% colnames(strand),
    "start" %in% colnames(strand),
    "end" %in% colnames(strand),
    "class" %in% colnames(strand)
  ) %>% invisible()
  strand[, class := paste("State:", class)]
  strand[, sample_cell := paste(sample, "-", cell)]

  set_diff <- setdiff(unique(strand$sample_cell), unique(counts$sample_cell))
  if (length(set_diff) > 0) message("[Warning] Strand states and Counts differ in cells: ", set_diff)

  strand <- strand[chrom == CHROM]
}

### Check complex regions file
if (!is.null(f_complex)) {
  message(" * Reading complex regions state file from ", f_complex, "...")
  complex <- fread(f_complex)
  assert_that(
    "chrom" %in% colnames(complex),
    "start" %in% colnames(complex),
    "end" %in% colnames(complex)
  ) %>% invisible()

  complex <- complex[chrom == CHROM]
  message("   --> Found ", nrow(complex), " complex region(s) in chromosome ", CHROM)
}

### Check SV groups file
if (!is.null(f_groups)) {
  message(" * Reading SV group file from ", f_groups, "...")
  groups <- fread(f_groups)
  assert_that(
    "chrom" %in% colnames(groups),
    "start" %in% colnames(groups),
    "end" %in% colnames(groups),
    "group_id" %in% colnames(groups)
  ) %>% invisible()
  groups[, group_id := paste("SV group", group_id)]
  groups <- groups[chrom == CHROM]
  message("   --> Found ", nrow(groups), " SV groups in chromosome ", CHROM)
}

### Check single cell segmentation file
if (!is.null(f_scsegments)) {
  message(" * Reading per-cell segmentation regions state file from ", f_scsegments, "...")
  scsegments <- fread(f_scsegments)
  assert_that(
    "sample" %in% colnames(scsegments),
    "cell" %in% colnames(scsegments),
    "chrom" %in% colnames(scsegments),
    "position" %in% colnames(scsegments)
  ) %>% invisible()
  scsegments[, sample_cell := paste(sample, "-", cell)]

  scsegments <- scsegments[chrom == CHROM]
}


################################################################################
# Actual plot                                                                  #
################################################################################


# Plot always a few cells per page!
y_lim <- 3 * counts[, median(w + c)]
n_cells <- length(unique(counts[, sample_cell]))
i <- 1


message(" * Plotting ", CHROM, " (", f_out, ")")
cairo_pdf(f_out, width = 14, height = 10, onefile = T)
while (i <= n_cells) {
  message(" * Processing cells from ", i, " to ", i + cells_per_page - 1)

  # Subset to this set of cells:
  CELLS <- unique(counts[, .(sample_cell)])[i:(min(i + cells_per_page - 1, n_cells))]
  setkey(CELLS, sample_cell)

  # Subset counts
  local_counts <- counts[CELLS, on = .(sample_cell), nomatch = 0]

  # Start major plot
  plt <- ggplot(local_counts)

  # Add background colors for segments, if available:
  if (!is.null(f_segments)) {
    message("   * Adding segment colors")
    # Segments need to be multiplied by "CELLS"
    local_seg <- CELLS[, as.data.table(seg), by = sample_cell]
    if (nrow(local_seg) > 0) {
      plt <- plt +
        geom_rect(
          data = local_seg, alpha = 0.4,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)
        )
    }
  }

  # Add colors for SV calls, if available
  if (!is.null(f_calls)) {
    message("   * Adding SV calls")
    local_svs <- svs[CELLS, on = .(sample_cell), nomatch = 0]
    if (nrow(local_svs) > 0) {
      plt <- plt +
        geom_rect(
          data = local_svs, alpha = 1,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = SV_class)
        )
    }
  }

  # Add bars for true SVs, if available
  if (!is.null(f_truth)) {
    message("   * Adding true SVs")
    local_sim <- simul[CELLS, on = .(sample_cell), nomatch = 0]
    if (nrow(local_sim) > 0) {
      plt <- plt +
        geom_rect(
          data = local_sim,
          aes(xmin = start, xmax = end, ymin = y_lim, ymax = Inf, fill = SV_class)
        )
    }
  }

  # Add lines for single cell segmentation, if available
  if (!is.null(f_scsegments)) {
    message("   * Adding single cell segments")
    local_scsegments <- scsegments[CELLS, on = .(sample_cell), nomatch = 0]
    if (nrow(local_scsegments) > 0) {
      plt <- plt +
        geom_segment(
          data = local_scsegments,
          aes(x = position, xend = position, y = -Inf, yend = -.8 * y_lim), color = "blue"
        )
    }
  }

  # Add bars for strand states, if available
  if (!is.null(f_strand)) {
    message("   * Adding strand states")
    local_strand <- strand[CELLS, on = .(sample_cell), nomatch = 0]
    if (nrow(local_strand) > 0) {
      plt <- plt +
        geom_rect(
          data = local_strand,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = -y_lim, fill = class)
        )
    }
  }

  # Add bars for SV group, if available
  if (!is.null(f_groups)) {
    message("   * Adding SV groups")
    if (nrow(groups) > 0) {
      plt <- plt +
        geom_rect(
          data = groups,
          aes(xmin = start, xmax = end, ymin = .85 * y_lim, ymax = Inf, fill = group_id)
        )
      # Add colors for SV classes
      manual_colors <- c(manual_colors, setNames(colorRampPalette(brewer.pal(12, "Set2"))(nrow(groups)), groups$group_id))
    }
  }

  # Add bars for complex states, if available
  if (!is.null(f_complex)) {
    message("   * Adding complex intervals")
    if (nrow(complex) > 0) {
      plt <- plt +
        geom_rect(
          data = complex,
          aes(xmin = start, xmax = end, ymin = y_lim, ymax = Inf), fill = "darkorchid1"
        )
    }
  }

  message("   * Adding actual W/C counts")
  plt <- plt +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = -w), fill = "sandybrown") +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = c), fill = "paleturquoise4")


  # Highlight None bins, if requested
  none_bins <- local_counts[class == "None"]
  if (show_none == T && nrow(none_bins) > 0) {
    message("   * Highlighting None bins")
    plt <- plt +
      geom_segment(data = none_bins, aes(x = start, xend = end, y = 0, yend = 0), col = "black", size = 2)
  }


  message("   * Adding labels, etc.")
  plt <- plt +
    facet_wrap(~sample_cell, ncol = 1) +
    ylab("Watson | Crick") + xlab(NULL) +
    scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
    scale_y_continuous(breaks = pretty_breaks(3)) +
    coord_cartesian(ylim = c(-y_lim, y_lim)) +
    scale_fill_manual(values = manual_colors) +
    theme_minimal() +
    theme(
      panel.spacing = unit(0, "lines"),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(color = "#eeeeee", fill = "#eeeeee"),
      strip.text = element_text(size = 5),
      legend.position = "bottom"
    ) +
    ggtitle(paste("data:", basename(f_counts), "chromosome:", CHROM))

  message("   * outputting")
  print(plt)
  i <- i + cells_per_page
} # while
dev.off()
