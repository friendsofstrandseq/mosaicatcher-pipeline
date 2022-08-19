suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))


args = commandArgs(trailingOnly = T)
if (length(args)!=3) {
  print("Usage: Rscript scale.R <count table> <norm factors> <out>")
  print("")
  print("       Normalize Strand-seq read counts. Divide the counts of all bins")
  print("       by a scaling factor (norm$scalar) and further black-list bins")
  print("       if requested in the normalizatio file (norm$class).")
  options( show.error.messages = F)
  stop()
}

# Read counts
message(" * Reading counts from ", args[1])
counts = fread(paste("zcat",args[1]))
assert_that(is.data.table(counts),
            "chrom"  %in% colnames(counts),
            "start"  %in% colnames(counts),
            "end"    %in% colnames(counts),
            "class"  %in% colnames(counts),
            "sample" %in% colnames(counts),
            "cell"   %in% colnames(counts),
            "w"      %in% colnames(counts),
            "c"      %in% colnames(counts)) %>% invisible
setkey(counts, chrom, start, end)

# Check that all cells have the same bins
bins <- unique(counts[, .(chrom, start, end)])
counts[,
       assert_that(all(.SD == bins), msg = "Not the same bins in all cells"),
       by = .(sample, cell),
       .SDcols = c("chrom", "start", "end")] %>% invisible


# remove bad cells
bad_cells <- counts[class == "None", .N, by = .(sample, cell)][N == nrow(bins)]
if (nrow(bad_cells)>0) {
  message(" * Removing ", nrow(bad_cells), " cells because thery were black-listed.")
  counts <- counts[!bad_cells, on = c("sample","cell")]
}

# Check that the "None" bins are all the same across cells
none_bins <- unique(counts[!bad_cells, on = c("sample","cell")][class == "None", .(chrom, start, end)])
if (nrow(none_bins) > 0) {
  counts[!bad_cells, on = c("sample","cell")][class == "None",
         assert_that(all(.SD == none_bins, msg = "None bins are not the same in all cells (excl. bad cells)")),
         by = .(sample, cell),
         .SDcols = c("chrom", "start", "end")] %>% invisible
}


# Read normalization factors
message(" * Reading norm file from ", args[2])
norm = fread(args[2])
assert_that(is.data.table(norm),
            "chrom"  %in% colnames(norm),
            "start"  %in% colnames(norm),
            "end"    %in% colnames(norm),
            "scalar" %in% colnames(norm)) %>% invisible
if ("class" %in% colnames(norm)) {
  norm <- norm[, .(chrom, start, end, scalar, norm_class = class)]
} else {
  norm <- norm[, .(chrom, start, end, scalar, norm_class = "good")]
}
setkey(norm, chrom, start, end)


# Set particular values of the norm_class to "None":
norm[scalar < 0.01, norm_class := "None"]


# annotate counts with scaling factor
counts <- merge(counts,
                norm,
                by = c("chrom","start","end"),
                all.x = T)

if (any(is.na(counts$scalar))) {
  message(" * Assign scalars: Could not match ",
          unique(counts[,.(chrom, start, end, scalar)])[is.na(scalar), .N],
          " bins (out of ",
          unique(counts[,.(chrom, start, end)])[,.N],
          ") -> set those to 1")
}

# Fill gaps in the norm file
counts[is.na(scalar), `:=`(scalar = 1, norm_class = "good")]

# Black-listing bins
test <- counts[!bad_cells, on = c("sample","cell")][cell == unique(cell)[1]]
test <- test[, .(count_None = sum(class == "None"),
         norm_None  = sum(norm_class == "None"),
         final_None = sum(class == "None" | norm_class == "None"))]
message(" * ", test$count_None, " bins were already black-listed; ", test$norm_None, " are blacklisted via the normalization, leading to a total of ", test$final_None)
counts[norm_class == "None", class := "None"]



# Apply normalization factor
counts[, `:=`(c = as.numeric(c), w = as.numeric(w))]
counts[class != "None", `:=`(c = c * scalar,
                             w = w * scalar)]
counts[class == "None", `:=`(c = 0.0, w = 0.0)]

message(" * Applying normalization: min = ",
        round(min(counts[class!="None", scalar]),3),
        ", max = ",
        round(max(counts[class!="None", scalar]), 3),
        ", median = ",
        median(unique(counts[,.(chrom,start,end,class,scalar)][class!="None", scalar])))


# Remove column
counts[, norm_class := NULL]
counts[, scalar := NULL]


# Write down table
message(" * Write data to ", args[3])
gz1 <- gzfile(args[3], "w")
write.table(counts, gz1, sep = "\t", quote = F, col.names = T, row.names =F)
close(gz1)