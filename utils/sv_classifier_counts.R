library(data.table)
library(assertthat)


# add_seg_counts
# Given segments + a large count table, calculate W/C counts per segment and cell.
#
# df = table with segments / cells.
# counts = table with raw counts
#
# returns an updated df
#
# Internally, create a [bins x cells] matrix to access cumulative counts quickly
add_seg_counts <- function(df, count_tab) {
  
    assert_that(is.data.table(df))
    assert_that("sample" %in% colnames(df),
                "cell"   %in% colnames(df),
                "chrom"  %in% colnames(df),
                "from"   %in% colnames(df))
  
    counts <- count_tab # copy
    assert_that("chrom" %in% colnames(counts),
                "start" %in% colnames(counts),
                "end"   %in% colnames(counts),
                "sample"%in% colnames(counts),
                "cell"  %in% colnames(counts))
    counts <- counts[order(sample,cell,chrom,start,end),] # order
    setkey(counts,sample,cell)
    
    

    # Assign bin indices and check that all cells have the same bins!
    counts[, idx := 1:.N, by = .(chrom, sample, cell)]
    bins <- unique(counts[, .(chrom, start, end)])
    counts[, assert_that(all(.(chrom,start,end) == bins)), by = .(sample,cell)]

    # chrom_map (starts with 0)
    chrom_map <- bins[,.N, by = chrom][,.(chrom = c(chrom, "end"), chr_idx = 1:(length(chrom)+1), N = c(1,cumsum(N)+1))]
    setkey(chrom_map, chrom)

    
    # Get cumulative counts (bins are ordered by idx!)
    counts[, w := cumsum(w), by = .(sample,cell,chrom)]
    counts[, c := cumsum(c), by = .(sample,cell,chrom)]

    # Make cumulative count MATRIX
    cC = dcast(counts, chrom + idx ~ sample + cell, value.var = "c")
    cW = dcast(counts, chrom + idx ~ sample + cell, value.var = "w")
    cC[, assert_that(!is.unsorted(idx)), by = chrom]
    cW[, assert_that(!is.unsorted(idx)), by = chrom]


    # Function to quickly access counts in any segemnt [from_, to_] 
    # (boundaries inclusive). Note that bin start with 1 in R.
    # accesses variables chrom_map, cW, and cC
    count <- function(sample_, cell_, chrom_, from_, to_) {
        assert_that(from_ > 0, 
                    to_ > 0,
                    from_ <= to_,
                    chrom_ %in% chrom_map$chrom,
                    to_ <= chrom_map$N[chrom_map[chrom_, chr_idx]+1])
        col.name = paste(sample_, cell_, sep = "_")
        assert_that(col.name %in% colnames(cW))
        offset = chrom_map[chrom_]$N -1
        assert_that(cW$idx[offset + to_] == to_)
        if (from_ == 1) {
          watson = cW[[col.name]] [offset + to_]
          crick  = cC[[col.name]] [offset + to_]
        } else {
          watson = cW[[col.name]] [offset + to_] - cW[[col.name]][offset + from_-1]
          crick  = cC[[col.name]] [offset + to_] - cC[[col.name]][offset + from_-1]
        }
        return (list(watson,crick))
    }

    
    df[, c("W","C") := count(sample, cell, chrom, from, to), by = .(sample, cell, chrom, from)  ]
    return (df)
}
