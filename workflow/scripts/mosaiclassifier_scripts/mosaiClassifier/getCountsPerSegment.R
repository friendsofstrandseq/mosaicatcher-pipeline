suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))


# addCountsPerSegment
# Given segments + a large count table, calculate W/C counts per segment and cell.
#
# df = table with segments / cells.
# counts = table with raw counts
#
# returns an updated df
#
# Internally, create a [bins x cells] matrix to access cumulative counts quickly
#
# Todo:
# A maybe faster way of doing it would be a more vectorized / data.table conform
# way. This would be a nice step for later
#
addCountsPerSegment2 <- function(df, count_tab) {
  assert_that(is.data.table(df))
  assert_that(
    "sample" %in% colnames(df),
    "cell" %in% colnames(df),
    "chrom" %in% colnames(df),
    "from" %in% colnames(df)
  )

  counts <- copy(count_tab) # copy
  assert_that(
    "chrom" %in% colnames(counts),
    "start" %in% colnames(counts),
    "end" %in% colnames(counts),
    "sample" %in% colnames(counts),
    "cell" %in% colnames(counts),
    "class" %in% colnames(counts)
  )
  counts <- counts[order(sample, cell, chrom, start, end), ] # order
  setkey(counts, sample, cell, chrom, start, end)

  # Add expected counts (old way)
  df[,
    expected_old := (to - from + 1) * mean,
    by = .(sample, cell, chrom, from, to)
  ]

  # Add expected counts (new way)
  counts[, num_bins := cumsum(class != "None"), by = .(sample, cell, chrom)]
  xxx <- counts[, .(chrom_ = chrom, start_ = start, end_ = end, sample_ = sample, cell_ = cell, num_bins)]
  df[,
    expected := xxx[sample_ == sample & cell_ == cell & chrom_ == chrom, num_bins[to] - num_bins[from] + 1] * mean,
    by = .(sample, cell, chrom)
  ]

  # unique(probs[, .(chrom, from, to)])

  # Assign bin indices and check that all cells have the same bins!
  counts[, idx := 1:.N, by = .(chrom, sample, cell)]
  bins <- unique(counts[, .(chrom, start, end)])
  counts[, assert_that(all(.(chrom, start, end) == bins)), by = .(sample, cell)]

  # chrom_map (starts with 0)
  chrom_map <- bins[, .N, by = chrom][, .(chrom = c(chrom, "end"), chr_idx = 1:(length(chrom) + 1), N = c(1, cumsum(N) + 1))]
  setkey(chrom_map, chrom)

  # Set black-listed counts to 0
  counts[class == "None", `:=`(w = 0, c = 0)]

  # Get cumulative counts (bins are ordered by idx!)
  counts[, w := cumsum(w), by = .(sample, cell, chrom)]
  counts[, c := cumsum(c), by = .(sample, cell, chrom)]

  # Make cumulative count MATRIX
  cC <- dcast(counts, chrom + idx ~ sample + cell, value.var = "c")
  cW <- dcast(counts, chrom + idx ~ sample + cell, value.var = "w")
  cC[, assert_that(!is.unsorted(idx)), by = chrom]
  cW[, assert_that(!is.unsorted(idx)), by = chrom]


  # Function to quickly access counts in any segemnt [from_, to_]
  # (boundaries inclusive). Note that bin start with 1 in R.
  # accesses variables chrom_map, cW, and cC
  count <- function(sample_, cell_, chrom_, from_, to_) {
    assert_that(
      from_ > 0,
      to_ > 0,
      from_ <= to_,
      chrom_ %in% chrom_map$chrom,
      to_ <= chrom_map$N[chrom_map[chrom_, chr_idx] + 1]
    )
    col.name <- paste(sample_, cell_, sep = "_")
    assert_that(col.name %in% colnames(cW))
    offset <- chrom_map[chrom_]$N - 1
    assert_that(cW$idx[offset + to_] == to_)
    if (from_ == 1) {
      watson <- cW[[col.name]][offset + to_]
      crick <- cC[[col.name]][offset + to_]
    } else {
      watson <- cW[[col.name]][offset + to_] - cW[[col.name]][offset + from_ - 1]
      crick <- cC[[col.name]][offset + to_] - cC[[col.name]][offset + from_ - 1]
    }
    return(list(watson, crick))
  }


  df[, c("W", "C") := count(sample, cell, chrom, from, to), by = .(sample, cell, chrom, from)]
  return(df)
}



addCountsPerSegment <- function(df, counts) {
  assert_that(is.data.table(df))
  assert_that(
    "sample" %in% colnames(df),
    "cell" %in% colnames(df),
    "chrom" %in% colnames(df),
    "from" %in% colnames(df)
  ) %>% invisible()

  count_tab <- copy(counts) # copy
  assert_that(
    "chrom" %in% colnames(count_tab),
    "start" %in% colnames(count_tab),
    "end" %in% colnames(count_tab),
    "sample" %in% colnames(count_tab),
    "cell" %in% colnames(count_tab),
    "class" %in% colnames(count_tab)
  ) %>% invisible()
  setkey(count_tab, sample, cell, chrom, start, end)

  # Set black-listed counts to 0
  count_tab[class == "None", `:=`(w = 0, c = 0)]

  # Assign bins (from count_tab) to segments
  # Extend segments by all numbers from `from` to `to`
  print(df)
  all_segs <- unique(df[, .(chrom, from, to)])[, .(bin = from:to), by = .(chrom, from, to)]
  print(all_segs)
  # Assign indices to bins in count_tab
  count_tab[, bin := 1:.N, by = .(sample, cell, chrom)]
  # Merge count_tab and segments by their bin IDs
  print(count_tab)
  count_tab <- merge(count_tab, all_segs, by = c("chrom", "bin"), all.x = T)
  print(count_tab)
  assert_that(all(!is.na(count_tab$from)), msg = "Segments should cover all bins") %>% invisible()

  # Now summarize counts and expectation per cell and segment
  count_tab <- count_tab[,
    .(
      C = as.integer(sum(c)),
      W = as.integer(sum(w)),
      num_bins = sum(class != "None")
    ),
    by = .(sample, cell, chrom, from, to)
  ]

  # Add information to the original `df` table
  df <- merge(df, count_tab, by = c("sample", "cell", "chrom", "from", "to"))
  df[, expected := num_bins * mean]

  return(df)
}




addNormalizationScalar <- function(df, counts, normVector) {
  assert_that(
    is.data.table(df),
    "chrom" %in% colnames(df),
    "start" %in% colnames(df),
    "end" %in% colnames(df),
    "sample" %in% colnames(df),
    "cell" %in% colnames(df)
  ) %>% invisible()
  setkey(df, sample, cell, chrom, start, end)
  assert_that(
    is.data.table(normVector),
    "chrom" %in% colnames(normVector),
    "start" %in% colnames(normVector),
    "end" %in% colnames(normVector),
    "scalar" %in% colnames(normVector)
  )
  assert_that(nrow(normVector) == nrow(normVector[, .(chrom, start, end)]))
  setkey(normVector, chrom, start, end)

  x <- merge(unique(counts[, .(chrom, start, end)]),
    normVector[, .(chrom, start, end, scalar)],
    by = c("chrom", "start", "end"),
    all.x = T
  )
  setkey(x, chrom, start, end)

  if (nrow(x[is.na(scalar)]) > 0) {
    message("[MosaiClassifier] WARNING: Normalization could not be set for ", nrow(x[is.na(scalar)]), " bins. Defaulting to 1")
  }

  getScalarSum <- function(chrom_, from_, to_) {
    return(x[chrom == chrom_, sum(scalar[from_:to_]) / (to_ - from_ + 1)])
  }
  df[,
    scalar := getScalarSum(chrom, from, to),
    by = .(chrom, from, to)
  ]
}