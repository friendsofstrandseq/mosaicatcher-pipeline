library(data.table)

MAX_NEGATIVE_VALUE = -1000


### REF
calc_p_ref_wc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p / (1-df$nb_p) / 2
  r_c = r_w
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_ref_ww = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_ref_cc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}

add_p_ref <- function(df, alpha = 0.05) {
  df$p_ref = MAX_NEGATIVE_VALUE
  df[state == "WW"]$p_ref = calc_p_ref_ww(df[state == "WW"], alpha)
  df[state == "WC"]$p_ref = calc_p_ref_wc(df[state == "WC"], alpha)
  df[state == "CW"]$p_ref = calc_p_ref_wc(df[state == "CW"], alpha)
  df[state == "CC"]$p_ref = calc_p_ref_cc(df[state == "CC"], alpha)
  df[p_ref < MAX_NEGATIVE_VALUE]$p_ref = MAX_NEGATIVE_VALUE
  df
}




### HOM_INV
calc_p_homInv_wc = function(df, alpha = 0.05) {
  return(calc_p_ref_wc(df, alpha))
}
calc_p_homInv_ww = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_homInv_cc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}

add_p_homInv <- function(df, alpha = 0.05) {
  df$p_homInv = MAX_NEGATIVE_VALUE
  df[state == "WW"]$p_homInv = calc_p_homInv_ww(df[state == "WW"], alpha)
  df[state == "WC"]$p_homInv = calc_p_homInv_wc(df[state == "WC"], alpha)
  df[state == "CW"]$p_homInv = calc_p_homInv_wc(df[state == "CW"], alpha)
  df[state == "CC"]$p_homInv = calc_p_homInv_cc(df[state == "CC"], alpha)
  df[p_homInv < MAX_NEGATIVE_VALUE]$p_homInv = MAX_NEGATIVE_VALUE
  df
}




### HET_INV
calc_p_hetInv_wc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  prob_h1 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  r_w = df$scalar * df$expected * (1-alpha) * df$nb_p /(1-df$nb_p)
  r_c = df$scalar * df$expected *   alpha   * df$nb_p /(1-df$nb_p)
  prob_h2 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  return(log( (exp(prob_h1) + exp(prob_h2))/2))
}
calc_p_hetInv_ww = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_hetInv_cc = function(df, alpha = 0.05) {
  return(calc_p_hetInv_ww(df,alpha))
}

add_p_hetInv <- function(df, alpha = 0.05) {
  df$p_hetInv = MAX_NEGATIVE_VALUE
  df[state == "WW"]$p_hetInv = calc_p_hetInv_ww(df[state == "WW"], alpha)
  df[state == "WC"]$p_hetInv = calc_p_hetInv_wc(df[state == "WC"], alpha)
  df[state == "CW"]$p_hetInv = calc_p_hetInv_wc(df[state == "CW"], alpha)
  df[state == "CC"]$p_hetInv = calc_p_hetInv_cc(df[state == "CC"], alpha)
  df[p_hetInv < MAX_NEGATIVE_VALUE]$p_hetInv = MAX_NEGATIVE_VALUE
  df
}




### HET_DEL
calc_p_hetDel_wc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  prob_h1 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  prob_h2 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  return(log((exp(prob_h1) + exp(prob_h2))/2))
}
calc_p_hetDel_ww = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_hetDel_cc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}

add_p_hetDel <- function(df, alpha = 0.05) {
  df$p_hetDel = MAX_NEGATIVE_VALUE
  df[state == "WW"]$p_hetDel = calc_p_hetDel_ww(df[state == "WW"], alpha)
  df[state == "WC"]$p_hetDel = calc_p_hetDel_wc(df[state == "WC"], alpha)
  df[state == "CW"]$p_hetDel = calc_p_hetDel_wc(df[state == "CW"], alpha)
  df[state == "CC"]$p_hetDel = calc_p_hetDel_cc(df[state == "CC"], alpha)
  df[p_hetDel < MAX_NEGATIVE_VALUE]$p_hetDel = MAX_NEGATIVE_VALUE
  df
}





### HOM DEL
calc_p_homDel_all = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}

add_p_homDel <- function(df, alpha = 0.05) {
  df$p_homDel = MAX_NEGATIVE_VALUE
  df[state %in% c("CC","WC","CW","WW")]$p_homDel = calc_p_homDel_all(df[state %in% c("CC","WC","CW","WW")], alpha)
  df
}






### HET_DUP
calc_p_hetDup_wc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * (1-alpha)
  prob_h1 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * (1-alpha)
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) / 2
  prob_h2 = log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p))
  return(log((exp(prob_h1) + exp(prob_h2))/2))
}
calc_p_hetDup_ww = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * 3 / 2 * (1-alpha)
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}
calc_p_hetDup_cc = function(df, alpha = 0.05) {
  r_w = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * alpha
  r_c = df$scalar * df$expected * df$nb_p /(1-df$nb_p) * 3 / 2 * (1-alpha)
  return(log(dnbinom(df$W, r_w, df$nb_p)) + log(dnbinom(df$C, r_c, df$nb_p)))
}

add_p_hetDup <- function(df, alpha = 0.05) {
  df$p_hetDup = MAX_NEGATIVE_VALUE
  df[state == "WW"]$p_hetDup = calc_p_hetDup_ww(df[state == "WW"], alpha)
  df[state == "WC"]$p_hetDup = calc_p_hetDup_wc(df[state == "WC"], alpha)
  df[state == "CW"]$p_hetDup = calc_p_hetDup_wc(df[state == "CW"], alpha)
  df[state == "CC"]$p_hetDup = calc_p_hetDup_cc(df[state == "CC"], alpha)
  df[p_hetDup < MAX_NEGATIVE_VALUE]$p_hetDup = MAX_NEGATIVE_VALUE
  df
}


### Summary function
add_NB_probs <- function(df, alpha = 0.05) {
  assert_that(is.data.table(df))
  assert_that("state"    %in% colnames(df),
              "W"        %in% colnames(df),
              "C"        %in% colnames(df),
              "scalar"   %in% colnames(df),
              "expected" %in% colnames(df),
              "nb_p"     %in% colnames(df))
  
  df <- add_p_ref(df, alpha)
  df <- add_p_homInv(df, alpha)
  df <- add_p_hetInv(df, alpha)
  df <- add_p_hetDel(df, alpha)
  df <- add_p_homDel(df, alpha)
  df <- add_p_hetDup(df, alpha)
  df
}