Pred_PLS_R <- function(xtrain1, ytrain, xtest1, lv) {
    source("workflow/scripts/scNOVA_scripts/script_PLSDA/pls_R_scNOVA.R")
    result_pls <- pls_R(xtrain1, ytrain, lv)
    B <- matrix(0, lv, lv)
    that1 <- matrix(0, 1, lv)
    for (l in 1:lv) {
        B[l, l] <- result_pls$pls_b[l]
        that1[1, l] <- t(as.matrix(xtest1)) %*% as.matrix(result_pls$pls_w[, l])
    }

    for (j in 1:lv) {
        ypred <- that1[, 1:j] %*% B[1:j, 1:j] %*% t(result_pls$pls_q[, 1:j])
    }
    return(ypred)
}
