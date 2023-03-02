library(pracma)
library(progress)
library(foreach)
library(doParallel)

# vip_null_R <- function(x, y, perm, lv_for_vip) {
#     m <- nrow(y)
#     n <- ncol(y)
#     tmp1 <- matrix(0, ncol(x), perm)

#     for (i in 1:perm) {
#         ind <- sample(m, m, replace = FALSE)
#         X1_r <- x[ind, ]
#         result_pls_rand <- pls_R(X1_r, y, lv_for_vip)

#         w <- result_pls_rand$pls_w[, 1:lv_for_vip]
#         r2 <- result_pls_rand$pls_ssq[1:lv_for_vip, 4]
#         result_vip_rand <- vip_R_v2(w, r2)
#         tmp1[, i] <- result_vip_rand
#         # cat(paste0(i, ' '))
#     }
#     tmp1 <- Reshape(tmp1, ncol(x) * perm, 1)
#     return(tmp1)
# }



vip_null_R <- function(x, y, perm, lv_for_vip) {
    m <- nrow(y)
    n <- ncol(y)
    tmp1 <- matrix(0, ncol(x), perm)

    opts <- list(progress = "text", verbose = FALSE)

    cl <- makeCluster(64, outfile = "cluster_vip_null_R.log")
    registerDoParallel(cl)

    foreach(i = 1:perm, .combine = "cbind", .packages = c("pracma"), .options.snow = opts) %dopar% {
        source("workflow/scripts/scNOVA_scripts/script_PLSDA/pls_R_scNOVA.R")
        source("workflow/scripts/scNOVA_scripts/script_PLSDA/vip_R_v2_scNOVA.R")
        ind <- sample(m, m, replace = FALSE)
        X1_r <- x[ind, ]
        result_pls_rand <- pls_R(X1_r, y, lv_for_vip)

        w <- result_pls_rand$pls_w[, 1:lv_for_vip]
        r2 <- result_pls_rand$pls_ssq[1:lv_for_vip, 4]
        result_vip_rand <- vip_R_v2(w, r2)
        result_vip_rand
    } -> tmp1

    # Stop the parallel backend
    stopCluster(cl)

    tmp1 <- Reshape(tmp1, ncol(x) * perm, 1)
    return(tmp1)
}
