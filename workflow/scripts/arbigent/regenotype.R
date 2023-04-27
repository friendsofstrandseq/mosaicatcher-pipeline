# Whoeps, 31 May 2020
# This is a perfectly input-output-oriented script.

# Input: - probabilities.R
#        - bed file with groups
#        - outdir
#        - CN_mapabilty per segment
#
# Output: - table with classifications, respective CN and mapability
#         - dumbbellplot
#         - beeswarmplots
# setwd('/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/scripts/pipeline/utils/regenotyper')
print("Initialising Regenotyper...")

# supposed to suppress warnings, but not working
oldw <- getOption("warn")
options(warn = -1)


# !/usr/bin/env Rscript
suppressMessages(library("optparse"))
suppressMessages(library("tidyr"))
suppressMessages(library("stringr"))


# INPUT INSTRUCTIONS
option_list <- list(
  make_option(c("-f", "--file"),
    type = "character", default = NULL,
    help = "probabilities.R, produced by mosaicatcher", metavar = "character"
  ),
  make_option(c("-b", "--bed"),
    type = "character", default = NULL,
    help = "a bed file specifiying labels/groups for the segments.", metavar = "character"
  ),
  make_option(c("-o", "--outdir"),
    type = "character", default = "./outputcorr/",
    help = "output dir name [default= %default]", metavar = "character"
  ),
  make_option(c("-c", "--cn_map"),
    type = "character", default = NULL,
    help = "average copy numbers and mapability for all segments in given bed file", metavar = "character"
  ),
  make_option(c("-p", "--path_to_self"),
    type = "character", default = NULL,
    help = "Path to regenotyper, if run from another directory.", metavar = "character"
  ),
  make_option(c("-m", "--mode"),
    type = "character", default = "",
    help = "Can be bulk or single-cell.", metavar = "character"
  )
  #  make_option(c("-s", "--sample_sex"), type='character', default=NULL,
  #              help="Needed for proper normalization of read counts in gonosomes", metavar="character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
path_to_regenotyper <- opt$path_to_self
if (is.null(path_to_regenotyper) == FALSE) {
  print(paste0("Going to path ", path_to_regenotyper))
  setwd(path_to_regenotyper)
} else {
  print(paste0("Staying in path", getwd()))
}

suppressMessages(source("workflow/scripts/arbigent/probability_helpers.R"))
suppressMessages(source("workflow/scripts/arbigent/regenotype_helpers.R"))
suppressMessages(source("workflow/scripts/arbigent/bulk_helpers.R"))








# Ok lets go! Calm the minds of impatient humans first of all.
print("Processing and summarizing information, making plots")
print("This can take a few minutes.")

p_link <- opt$file
labels_link <- NULL # opt$bed
outdir_raw <- opt$outdir
debug_file <- opt$cn_map
# sample_sex = opt$sample_sex
m <- opt$mode

if (m == "single-cell") {
  ### switch modules on/off
  make_bell_bulk <- F
  make_table_bulk <- F
  make_bee_bulk <- F

  make_bell_sc <- F
  make_table_sc <- F
  make_bee_sc <- F

  run_singlecell_mode <- T
} else {
  make_bell_bulk <- T
  make_table_bulk <- T
  make_bee_bulk <- F

  make_bell_sc <- F
  make_table_sc <- F
  make_bee_sc <- F

  run_singlecell_mode <- F
}
suppressMessages(dir.create(outdir_raw))

# is input file specified?
if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("Please specify path to probabilities.R", call. = FALSE)
}

# Tell the user if we go for bed mode or single mode. At the occasion, also load it.
if (is.null(labels_link)) {
  print("No bed file specified. Examining everything together.")
  labels <- NULL
} else {
  print("Bed file provided. Will use it to stratify results w.r.t. groups")
  labels <- read.table(labels_link, stringsAsFactors = FALSE)
  colnames(labels) <- c("chrom", "start", "end", "group")

  # Bed files can be weird. For safety, we remove duplicate lines
  labels <- unique(labels)
}

# sample name could be interesting
# print(p_link)
sname <- tail(strsplit(sub("/arbigent_mosaiclassifier/sv_probabilities/probabilities.Rdata", "", p_link), "/")[[1]], 1)
# sname <- str_match(p_link, "([HG|NA|GM]+[0-9]{3,5})")[, 2]
# stop()

# load p to p_grouped
print("Loading probabilities table")
probs_raw <- load_and_prep_pdf(p_link)
print("Done loading")

print(unique(probs_raw$chrom))
# load file containing info on valid bins
CN <- read.table(debug_file, header = 1, stringsAsFactors = F)

# Cut CN down to important cols, then join it with probs so we know which inversion has how
# many valid bins
if (!is.null(CN)) {
  CNmerge <- as.data.frame(lapply(CN[, c("chrom", "start", "end", "valid_bins")], as.character))

  CNmerge <- as.tbl(CNmerge)
  CNmerge <- CNmerge %>% mutate(
    chrom = as.character(chrom),
    start = as.numeric(as.character(start)),
    end = as.numeric(as.character(end))
  ) # ,
  p2 <- full_join(probs_raw, CNmerge, by = c("chrom", "start", "end"))
}

# Remove invalid bins from probs_raw
probs_raw <- p2[!is.na(p2$sample), ]
# probs_raw = probs_raw[probs_raw$chrom == 'chr22',]

print(unique(probs_raw$chrom))
#############################################################
#############################################################
# ATTENTION LADIES AND GENTLEMEN, HERE IS THE NORMALIZATION #
# PLEASE PAY CLOSE ATTENTION TO THIS                        #
#############################################################

# Using the valid_bin info, we can reconstruct the length
# normalization factor here. It is between 0 and 1.
# Basically it is the mapability ratio (0 nothing, 1 perfect)
len_normalization <- as.numeric(as.character(probs_raw$valid_bins)) /
  ((probs_raw$end - probs_raw$start) / 100.)
len_normalization[len_normalization == 0] <- 1

# Depending on how manual segment counts were normalized before,
# and depending on which normalization you want to have, you
# have to choose different options here.
# If they have been length normalized, I am using both options
# A and B to arrive at the 'downscaling' solution ('2'). If A and
# B are disabled, no further normalization is done. Since W and C
# counts have been up-scaled by watson_crick_counts.py, this is
# also ok ('solution 1')
# By all means check the qc plots that arrive at the end of the
# snakemake!

# option A: if we multiply by len_normalization, we remove previous
# normalization, and return back to non-length-corrected counts.
# probs_raw$W <- probs_raw$W * len_normalization
# probs_raw$C <- probs_raw$C * len_normalization

# option B: we can downscale expectations
probs_raw$expected <- probs_raw$expected * len_normalization


##############################################################
##############################################################
##############################################################
# Additionally: adjust expectations based on biological sex ##
##############################################################

# if (sample_sex == 'male'){

# UPDATE: NO WE DONT HAVE TO NORMALIZE.
# In male samples, we expect half the number of reads on X ...
# probs_raw[probs_raw$chrom=='chrX',]$expected = (probs_raw[probs_raw$chrom=='chrX',]$expected) / 2
# And half in y
# probs_raw[probs_raw$chrom=='chrY',]$expected = (probs_raw[probs_raw$chrom=='chrY',]$expected) / 2

# Remove WC and CW cells in chrX and Y.
#  probs_raw = probs_raw[!((probs_raw$chrom == 'chrX') & (probs_raw$class %in% c('WC','CW'))),]
#  probs_raw = probs_raw[!((probs_raw$chrom == 'chrY') & (probs_raw$class %in% c('WC','CW'))),]
# } else {
#  probs_raw = probs_raw[!(probs_raw$chrom == 'chrY'),]
# }



# Adding group information to probs_raw.
if (is.null(labels)) {
  # If bed file was not provided, everyone is group 'all'
  probs_raw$group <- "all"
} else {
  # Else, the ones with a group get that one from the bed file
  probs_raw <- full_join(probs_raw, labels)
  # ... the remaining ones are called ungrouped.
  probs_raw$group[is.na(probs_raw$group)] <- "ungrouped"
}

# remove inf things [DIRTY SOLUTION! SHOULD BE DONE BETTER! Not needed apparently? Not sure. CHECK]
# probs_raw = probs_raw[probs_raw$logllh != 'Inf',]
# probs_raw = probs_raw[probs_raw$logllh != '-Inf',]
# probs_raw = probs_raw[!(probs_raw$cell %in% unique(probs_raw[probs_raw$logllh == '-Inf',]$cell)),]

###################################
##### OKAAYYY HERE WE GOOOOO ######
#### THIS IS THE MAIN WORKHORSE ###
###################################
print(unique(probs_raw$group))
# group = unique(probs_raw$group)[1] #for quick manual mode
for (group in unique(probs_raw$group)) {
  # Talk to human
  print(paste0("Running samples with group ", group))

  # Make the outfolder (maybe not necessary?)
  outdir <- gsub("\\.:", "_:", paste0(outdir_raw, group, "/"))
  suppressMessages(dir.create(outdir))

  # cut pg down to the desired inversions
  pg <- probs_raw[probs_raw$group == group, ]
  haps_to_consider <- na.omit(unique(pg$haplotype))


  if (make_bell_bulk) {
    ### [I]a) make dumbbell plot ###
    #### [I] BULK ###


    pg_bulk_list <- (bulkify_pg(haps_to_consider, pg))
    pg_bulk <- data.frame(pg_bulk_list[1]) %>% group_by(start, end, haplotype, class, group)
    pg_bulk_probs <- data.frame(pg_bulk_list[2]) %>% group_by(start, end, haplotype, class, group)
    # write.table(pg_bulk, file=paste0('/home/hoeps/Desktop/', 'counts_bulk.txt'), quote = F, row.names = F, col.names = T)

    # at least temporarily, I'm operating both with likelihoods and probabilities. Haven't decided yet which one I like more.
    # call_llhs_bulk = (make_condensed_sumlist(haps_to_consider, pg_bulk)) %>% mutate_all(funs(replace_na(.,-1000)))
    call_llhs_bulk <- (make_condensed_sumlist(haps_to_consider, pg_bulk)) %>% mutate_all(funs(replace_na(., -1000)))

    # call_probs_bulk = (make_condensed_sumlist_probs(haps_to_consider, pg_bulk_probs)) %>% mutate_all(funs(replace_na(.,-1000)))
    # write.table(mm, file=paste0('/home/hoeps/Desktop/', 'counts_bulk2.txt'), quote = F, row.names = F, col.names = T)



    g <- make_dumbbell(call_llhs_bulk, groupname = group, run_shiny = F)
    # p = suppressMessages(ggplotly(g))

    savepath <- paste0(outdir, "bellplot_bulk.html")
    # suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
    ggsave(filename = paste0(outdir, sname, "_", group, "_bellplot_bulk.png"), width = 30, height = 12, units = "cm", device = "png")
    ggsave(filename = paste0(outdir, sname, "_", group, "_bellplot_bulk.pdf"), width = 30, height = 12, units = "cm", device = "pdf")

    # call_probs_bulk[,4:73] = -log(1-(call_probs_bulk[,4:73]))
    # call_probs_bulk[,4:73] = 10**(call_probs_bulk[,4:73])
    # g = make_dumbbell_probs(call_probs_bulk, groupname=group, run_shiny=F)
    # p = suppressMessages(ggplotly(g))

    savepath <- paste0(outdir, "bellplot_bulk_prob.html")
    # suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)),basename(savepath))))
    ggsave(g, filename = paste0(outdir, sname, "_", group, "_bellplot_bulk_prob.png"), width = 30, height = 12, units = "cm", device = "png")
    ggsave(g, filename = paste0(outdir, sname, "_", group, "_bellplot_bulk_prob.pdf"), width = 30, height = 12, units = "cm", device = "pdf")
  }

  if (make_table_bulk) {
    #### [I]b) write table ####
    tab <- make_table_finaledition(call_llhs_bulk, group, sname)
    # write.table(tab, file=paste0('/home/hoeps/Desktop/', 'counts_bulk_labels.txt'), quote = F, row.names = F, col.names = T)
    # adding copy number and mapability information to the table
    tab2 <- left_join(tab, CN[, c("chrom", "start", "end", "valid_bins")])


    if (dim(tab2[tab2$valid_bins == 0, ])[1] > 0) {
      tab2[tab2$valid_bins == 0, ]$pred_hard <- "nomappability"
      tab2[tab2$valid_bins == 0, ]$pred_soft <- "nomappability"
      tab2[tab2$valid_bins == 0, ]$pred_nobias <- "nomappability"
      tab2[tab2$valid_bins == 0, ]$second_hard <- "nomappability"
      tab2[tab2$valid_bins == 0, ]$confidence_hard_over_second <- 100
    }
    tab <- tab2[, !(names(tab2) %in% c("valid_bins"))]
    # t2 = as.data.frame(lapply(tab, as.character))

    tab[tab$confidence_hard_over_second == 0, "pred_hard"] <- "0|0"


    # if (!is.null(CN)){
    # CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","CN","mapability")], as.character))
    #  CNmerge = as.data.frame(lapply(CN[, c("chrom","start","end","valid_bins")], as.character))

    #  tab <- left_join(t2, CNmerge, by = c("chrom","start","end"))
    # CN = read.table(CN_link, stringsAsFactors = F, header=1);
    # }

    tab$group <- group
    write.table(tab, file = paste0(outdir, "sv_calls_bulk.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  }

  if (make_bee_bulk) {
    ## [I]c) make beeswarm plots ##
    save_beeswarms(pg_bulk %>% group_by(start), call_llhs_bulk, outdir, testrun = F, compositemode = T)
  }

  ### [II] SINGLE CELL ###


  if (make_bell_sc) {
    call_llhs <- (make_condensed_sumlist(haps_to_consider, pg)) # %>% mutate_all(funs(replace_na(.,-1000)))

    #### [II]a) make dumbbell plot ####


    suppressMessages(source("workflow/scripts/arbigent/regenotype_helpers.R")) # for quick manual mode

    # create the ggplot plot
    g <- make_dumbbell(call_llhs, groupname = group, run_shiny = F)
    g
    # convert it to plotly
    # p = suppressMessages(ggplotly(g))

    # save. htmlwidgets does not work with relative paths, so we need a little workaround
    # (info and code taken from https://stackoverflow.com/questions/41399795/savewidget-from-htmlwidget-in-r-cannot-save-html-file-in-another-folder)
    savepath <- paste0(outdir, "bellplot.html")
    suppressMessages(htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(savepath)), basename(savepath))))
    ggsave(filename = paste0(outdir, sname, "_", group, "_bellplot.png"), width = 30, height = 12, units = "cm", device = "png")
    ggsave(filename = paste0(outdir, sname, "_", group, "_bellplot.pdf"), width = 30, height = 12, units = "cm", device = "pdf")
  }

  if (make_table_sc) {
    #### [II]b) write table ####
    # tab = make_table(call_llhs, group, sname)
    tab <- make_table_finaledition(call_llhs, group, sname)

    t2 <- as.data.frame(lapply(tab, as.character))
    CNmerge <- as.data.frame(lapply(CN[, c("chrom", "start", "end", "CN", "mapability")], as.character))
    tab2 <- left_join(t2, CNmerge, by = c("chrom", "start", "end"))


    # adding copy number and mapability information to the table
    write.table(tab2, file = paste0(outdir, "sv_calls.txt"), quote = F, row.names = F, col.names = T)
  }

  if (make_bee_sc) {
    #### [II]c) save beewarm plots ###
    save_beeswarms(pg, call_llhs, outdir, testrun = F)
    print(paste0("Group ", group, " done."))
  }


  if (run_singlecell_mode) {
    # rectify logllh Nan. It will be overwritten anyway. Just want to avoid NA errors
    if (dim(pg[is.na(pg$logllh), ])[1] > 0) {
      pg[is.na(pg$logllh), ]$logllh <- -1
    }

    # re-calculate logllh based on mapping parameters
    pg2 <- calc_new_logllhs_singlecell(pg)

    # make an output table

    tab <- make_table_sc_separated(pg2)
    cols_to_return <- c("cell", "chrom", "start", "end", "class", "expected", "W", "C", "top_pred", "second_pred", "llr_1st_to_2nd", "llr_1st_to_ref")
    tab <- make_table_sc_separated(pg2)
    write.table(tab[, cols_to_return], file = paste0(outdir_raw, "sv_calls.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
    write.table(tab, file = paste0(outdir_raw, "sv_calls_detailed.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  }
}

print("### ALL DONE. Happy discoveries. ###")

# Does not seem to work.
options(warn = oldw)



# Manual links for debugging.
# p_link = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/sv_probabilities/HG00733/100000_fixed_norm.selected_j0.1_s0.5/probabilities.Rdata'
# p_link ='../results/audano_calls/probabilities.Rdata'
# outdir = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/31st_may/'
