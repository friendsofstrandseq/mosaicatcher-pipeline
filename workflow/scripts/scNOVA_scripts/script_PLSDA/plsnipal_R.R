plsnipal_R <- function(x,y){
my <- nrow(y); ny <- ncol(y)
if (ny > 1){
  ssy <- round(colSums(y^2));
  ymax <- max(ssy); yi <- which.max(ssy)
  u <- y[,yi];
} else {
  u <- y[,1];
}
conv <- 1;
told <- x[,1];
count <- 1.0;

#  Specify the conversion tolerance
while (conv > 1e-4){
  count <- count + 1;
  w <- t(t(as.matrix(u)) %*% as.matrix(x));
  w <- t(t(w)/max(svd(t(w))$d));
  t <- as.matrix(x) %*% as.matrix(w);
  if (ny == 1){
    q = 1;
    break
  }
  q <- t(t(t) %*% as.matrix(y));
  q <- t(t(q)/max(svd(t(q))$d));
  u <- as.matrix(y) %*% as.matrix(q);
  conv = max(svd(told - t)$d);
  told = t;
  if (count >= 1000){
    cat('Algorithm Failed to Converge after 1000 Iterations')
    conv
    break
  }
}
tmp <- (as.matrix(t(t)) %*% as.matrix(t));
p <- t( (as.matrix(t(t)) %*% as.matrix(x)) / as.numeric(tmp) );
p_norm <- max(svd(p)$d);
t <- t*p_norm;
w = w*p_norm;
p = p/p_norm;

return(list(result_p = p, result_q = q, result_w = w, result_t = t, result_u = u))
}