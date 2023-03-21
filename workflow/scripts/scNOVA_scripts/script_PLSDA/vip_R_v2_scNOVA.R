vip_R_v2 <- function(w, r2) {
  k <- nrow(w)
  a <- ncol(w)
  # library(SuperPCA)
  source("workflow/scripts/scNOVA_scripts/script_PLSDA/normc.R")
  vip <- sqrt(rowSums((normc(w)^2) %*% diag(r2)) / sum(r2) * nrow(w))
  return(vip)
}
