suppressMessages(library(mc2d))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
source("workflow/scripts/mosaiclassifier_scripts/mosaiClassifier/getDispParAndSegType.R")

# defining multinomial parameters based on haplosegType and alpha
getMultinomialParams <- function(haploSegType, alpha) {
  CNs <- NULL
  if (haploSegType == "?") {
    CNs <- rep(-1, 4)
  } else {
    CNs <- as.numeric(unlist(strsplit(haploSegType, "")))
    CNs[CNs == 0] <- alpha
    CNs <- CNs / sum(CNs)
  }

  return(list(
    multi.p.W.h1 = CNs[1],
    multi.p.C.h1 = CNs[2],
    multi.p.W.h2 = CNs[3],
    multi.p.C.h2 = CNs[4]
  ))
}

# adding an extra column for multinomial probabilities
addHaploCountProbs <- function(probs, haploCounts, alpha) {
  assert_that(
    is.data.table(probs),
    "sample" %in% colnames(probs),
    "cell" %in% colnames(probs),
    "chrom" %in% colnames(probs),
    "start" %in% colnames(probs),
    "end" %in% colnames(probs),
    "nb_p" %in% colnames(probs),
    "expected" %in% colnames(probs),
    "W" %in% colnames(probs),
    "C" %in% colnames(probs),
    "class" %in% colnames(probs),
    "haplotype" %in% colnames(probs)
  ) %>% invisible()

  # sorting probs
  setkey(probs, sample, chrom, start, end, cell, haplotype)

  # adding segTypes and haplo segTypes
  probs[, `:=`(
    segtype = paste0(getSegType(class[1], haplotype[1]), collapse = ""),
    haploSegType = ifelse(class == "?", "?", paste0(c(
      getSegType(substr(class[1], 1, 1), substr(haplotype[1], 1, 2)),
      getSegType(substr(class[1], 2, 2), substr(haplotype[1], 3, 4))
    ), collapse = ""))
  ),
  by = .(class, haplotype)
  ]

  # merging probs and haplotype counts
  probs <- merge(probs,
    haploCounts[, .(chrom, start, end, cell, watson.H1, crick.H1, watson.H2, crick.H2)],
    by = c("chrom", "start", "end", "cell"),
    allow.cartesian = T,
    all.x = T
  )

  # define multinomial probabilities for each of the four different outcomes based on haplo segtype
  probs[, c("multi.p.W.h1", "multi.p.C.h1", "multi.p.W.h2", "multi.p.C.h2") := getMultinomialParams(haploSegType[1], alpha),
    by = haploSegType
  ]

  # add a column for multinomial probabilities of haplotagged read counts
  probs.new <- probs[!is.na(watson.H1) & !is.na(multi.p.W.h1) & haploSegType != "?"]

  probs.new[, haplotag.prob := dmultinomial(
    x = as.matrix(probs.new[, .(watson.H1, crick.H1, watson.H2, crick.H2)]),
    prob = as.matrix(probs.new[, .(multi.p.W.h1, multi.p.C.h1, multi.p.W.h2, multi.p.C.h2)])
  )]

  # merge probs.new and probs
  probs <- merge(probs,
    probs.new[, .(chrom, start, end, cell, sample, haplotype, haplotag.prob)],
    by = c("chrom", "start", "end", "cell", "sample", "haplotype"),
    all.x = T,
    allow.cartesian = T
  )

  return(probs)
}