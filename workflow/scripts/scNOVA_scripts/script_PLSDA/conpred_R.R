conpred_R <- function(b,w,p,q,lv){

mq <- nrow(q); nq <- ncol(q);
mw <- nrow(w); nw <- ncol(w);
if (nw != lv){
  if (lv > nw){
    cat(paste0('Original model has a maximum of ', nw,' LVs Calculating vectors for ', nw, ' LVs only'))
  	lv <- nw;
  } else {
    w <- w[,1:lv];
	  q <- q[,1:lv];
	  p <- p[,1:lv];
	  b <- b[,1:lv];
  }
}
m <- matrix(0, mq*lv, mw);
if (mq == 1){
  m <- t(as.matrix(w) %*% (solve(as.matrix(t(p)) %*% as.matrix(w))) %*% diag(as.vector(b)));
} else {
  mp <- t(as.matrix(w) %*% (solve(as.matrix(t(p)) %*% as.matrix(w))) %*% diag(as.vector(b)));
  for (i in 1:lv){
    mpp <- mp[i,];
    m[((i-1)*mq+1):(i*mq),] <- diag(q[,i]) %*% t(matrix(1, nrow(p), nrow(q))*mpp);
  }
}
return(m)
}