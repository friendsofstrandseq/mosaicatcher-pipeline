
# simplest GT assignments imaginable
add_gts <- function(tab, bias_factor, cutoff) {
  tab$GT <- "UNK"
  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GT <- tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias

  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias, "_lowconf")


  tab[tab$confidence_nobias_over_hard < bias_factor, ]$GT <- tab[tab$confidence_nobias_over_hard < bias_factor, ]$pred_hard

  tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    paste0(tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard, "_lowconf")

  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GT"] <- "noreads"
  return(tab)
}

add_gts_revisited <- function(tab, bias_factor, bias_add_factor, cutoff) {
  tab$GT <- "UNK"
  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GT <- tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias

  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias)


  tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$GT <-
    tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_hard

  tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard # removed lowconf here.

  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GT"] <- "noreads"
  return(tab)
}

add_gts_revisited_lowconf <- function(tab, bias_factor, bias_add_factor, cutoff) {
  tab$GT <- "UNK"
  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GT <- tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias
  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias)

  tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$GT <-
    tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_hard
  tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$GT <-
    paste0(tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard, "_lowconf") # removed lowconf here.

  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GT"] <- "noreads"
  return(tab)
}

add_long_gts_revisited_oldfmt <- function(tab, bias_factor, bias_add_factor, cutoff) {
  tab$GTL <- "UNK"

  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GTL <-
    paste(tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second + tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_nobias_over_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias),
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$GTL <-
    paste(tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_hard,
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_nobias,
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_hard_over_second,
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_nobias_over_hard,
      sep = ":"
    )

  tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(paste0(tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard, "_lowconf"),
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias,
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      sep = ":"
    )


  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GTL"] <- "noreads"
  return(tab)
}

add_long_gts_revisited <- function(tab, bias_factor, bias_add_factor, cutoff) {
  tab$GTL <- "UNK"

  tab$confidence_hard_over_second <- round(tab$confidence_hard_over_second, 1)
  tab$confidence_nobias_over_hard <- round(tab$confidence_nobias_over_hard, 1)


  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GTL <-
    paste(tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias,
      "T",
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second + tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_nobias_over_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias),
      "T",
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$GTL <-
    paste(tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_hard,
      "T",
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$pred_nobias,
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_hard_over_second,
      tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$confidence_nobias_over_hard,
      sep = ":"
    )

  tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard,
      "F",
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias,
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      sep = ":"
    )


  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GTL"] <- "noreads"
  return(tab)
}


add_long_gts <- function(tab, bias_factor, cutoff) {
  tab$GTL <- "UNK"

  tab[tab$confidence_nobias_over_hard >= bias_factor, ]$GTL <-
    paste(tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_nobias,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$pred_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second + tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_nobias_over_hard,
      tab[tab$confidence_nobias_over_hard >= bias_factor, ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(paste0(tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias, "_lowconf"),
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      sep = ":"
    )

  tab[tab$confidence_nobias_over_hard < bias_factor, ]$GTL <-
    paste(tab[tab$confidence_nobias_over_hard < bias_factor, ]$pred_hard,
      tab[tab$confidence_nobias_over_hard < bias_factor, ]$pred_nobias,
      tab[tab$confidence_nobias_over_hard < bias_factor, ]$confidence_hard_over_second,
      tab[tab$confidence_nobias_over_hard < bias_factor, ]$confidence_hard_over_second +
        tab[tab$confidence_nobias_over_hard < bias_factor, ]$confidence_nobias_over_hard,
      sep = ":"
    )

  tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GTL <-
    paste(paste0(tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_hard, "_lowconf"),
      tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$pred_nobias,
      tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second,
      tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_hard_over_second +
        tab[(tab$confidence_nobias_over_hard < bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$confidence_nobias_over_hard,
      sep = ":"
    )


  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GTL"] <- "noreads"
  return(tab)
}



add_simple_gts <- function(tab_f) {
  tab <- tab_f
  tab$GTs <- "UNK"
  tab[((tab$confidence_nobias_over_hard >= bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)), ]$GTs <- "complex"
  tab[(tab$confidence_nobias_over_hard >= bias_factor) & (tab$confidence_hard_over_second < cutoff), ]$GTs <- "complex"
  tab[(tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor), ]$GTs <- "simple"
  tab[((tab$confidence_nobias_over_hard < bias_factor) | (tab$confidence_nobias_over_hard < bias_add_factor)) & (tab$confidence_hard_over_second < cutoff), ]$GTs <- "simple_lowconf"
  tab[(tab$confidence_hard_over_second == 0) & (tab$confidence_nobias_over_hard == 0), "GTs"] <- "noreads"
  return(tab)
}
