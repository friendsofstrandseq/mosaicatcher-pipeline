pls_R <- function(x, y, lv) {
  out <- 1
  mx <- nrow(x)
  nx <- ncol(x)
  my <- nrow(y)
  ny <- ncol(y)
  if (nx < lv) {
    cat("No. of LVs must be <= no. of x-block variables")
  }
  p <- matrix(0, nx, lv)
  q <- matrix(0, ny, lv)
  w <- matrix(0, nx, lv)
  t <- matrix(0, mx, lv)
  u <- matrix(0, my, lv)
  b <- matrix(0, 1, lv)
  ssq <- matrix(0, lv, 2)
  ssqx <- sum(colSums(x^2))
  ssqy <- sum(colSums(y^2))
  olv <- lv
  rankx <- qr(x)$rank
  if (rankx < olv) {
    lv <- rankx
    if (out == 1) {
      cat(paste0("Rank of X is ", lv, ", which is less than lv of ", olv, " Calculating ", lv, " LVs only"))
    }
  }


  for (i in 1:lv) {
    # cat(paste0(i, ' '))
    source("workflow/scripts/scNOVA_scripts/script_PLSDA/plsnipal_R.R")
    result <- plsnipal_R(x, y)
    pp <- result$result_p
    qq <- result$result_q
    ww <- result$result_w
    tt <- result$result_t
    uu <- result$result_u

    b[1, i] <- (as.matrix(t(uu)) %*% as.matrix(tt)) / (as.matrix(t(tt)) %*% as.matrix(tt))
    x <- x - (as.matrix(tt) %*% as.matrix(t(pp)))
    y <- y - (b[1, i] * as.matrix(tt) %*% as.matrix(t(qq)))
    ssq[i, 1] <- (sum(colSums(x^2))) * 100 / ssqx
    ssq[i, 2] <- (sum(colSums(y^2))) * 100 / ssqy
    t[, i] <- tt[, 1]
    u[, i] <- uu[, 1]
    p[, i] <- pp[, 1]
    w[, i] <- ww[, 1]
    q[, i] <- qq[, 1]
  }

  if (olv > lv) {
    ssq[(lv + 1):nrow(ssq), 2] <- ssq[lv, 2]
  }
  ssqdif <- matrix(0, lv, 2)
  ssqdif[1, 1] <- 100 - ssq[1, 1]
  ssqdif[1, 2] <- 100 - ssq[1, 2]

  for (i in 2:olv) {
    for (j in 1:2) {
      ssqdif[i, j] <- -ssq[i, j] + ssq[i - 1, j]
    }
  }


  ssq <- cbind((1:olv), ssqdif[, 1], cumsum(ssqdif[, 1]), ssqdif[, 2], cumsum(ssqdif[, 2]))
  colnames(ssq) <- c("LV", "X_This_LV", "X_Total", "Y_This_LV", "Y_Total")


  m <- matrix(0, olv * ny, nx)
  source("workflow/scripts/scNOVA_scripts/script_PLSDA/conpred_R.R")
  m[1:(lv * ny), ] <- conpred_R(b, w, p, q, lv)


  if (ny > 1) {
    for (i in 2:olv) {
      j <- (i - 1) * ny + 1
      i0 <- j - ny
      m[j:(i * ny), ] <- m[j:(i * ny), ] + m[i0:((i - 1) * ny), ]
    }
  } else {
    for (k in 1:ncol(m)) {
      m[, k] <- cumsum(m[, k])
    }
  }


  return(list(pls_m = m, pls_ssq = ssq, pls_p = p, pls_q = q, pls_w = w, pls_t = t, pls_u = u, pls_b = b))
}
