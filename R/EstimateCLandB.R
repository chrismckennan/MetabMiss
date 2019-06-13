library(parallel)
######Estimator for C and L that is invariant to the missingness mechanism######


#####Estimate C from metabolites with <= t% missingness#####
EstC.0 <- function(Y, K, Cov=NULL, max.miss=0.05, max.iter=100, tol=1e-6, n.repeat.Sigma=3) {
  Frac.missing <- apply(Y, 1, function(y) {sum(is.na(y))/length(y)})
  ind.use <- Frac.missing <= max.miss
  Y <- Y[ind.use,]
  Frac.missing <- Frac.missing[ind.use]
  ind.no.miss <- which(Frac.missing == 0)
  ind.miss <- which(Frac.missing > 0)
  n <- ncol(Y)
  p <- nrow(Y)
  Ind.Obs <- !is.na(Y)
  
  if (!is.null(Cov)) {
    Cov <- cbind(Cov)
    Q <- Compute.Q(Cov)
    d <- ncol(Cov)
  } else {
    Q <- diag(n)
    d <- 0
  }
  C <- Q%*%cbind(svd(Y[Frac.missing==0,]%*%Q, nv=K)$v) * sqrt(n-d)
  W <- cbind(Cov, C)
  Sigma <- rep(1, p)
  BL <- matrix(NA, nrow=p, ncol=d+K)
  
  for (i in 1:n.repeat.Sigma) {
    if (i > 1) {
      Q.W <- Compute.Q(W)
      Sigma[Frac.missing==0] <- rowSums((Y[Frac.missing==0,]%*%Q.W)^2)/(n-ncol(W))
      if (sum(Frac.missing>0) > 0) {
        Sigma[Frac.missing > 0] <- unlist(lapply(ind.miss, function(g){ y<-Y[g,]; ind.obs<-!is.na(y); beta.hat <- solve(t(W[ind.obs,])%*%W[ind.obs,],t(W[ind.obs,])%*%y[ind.obs]); sum((y[ind.obs]-W[ind.obs,]%*%beta.hat)^2)/(sum(ind.obs)-ncol(W)) }))
      }
    }
    log.like.0 <- -1e16
    log.like.vec <- rep(NA, max.iter)
    
    for (j in 1:max.iter) {
      Q.W <- Compute.Q(W)
      
      ##Update BL##
      BL[Frac.missing==0,] <- Y[Frac.missing==0,]%*%W%*%solve(t(W)%*%W)
      if (sum(Frac.missing > 0) > 0) {
        BL[Frac.missing > 0,] <- t(apply(Y[Frac.missing > 0,], 1, function(y) {ind.obs <- !is.na(y); as.vector( solve(t(W[ind.obs,])%*%W[ind.obs,], t(W[ind.obs,])%*%y[ind.obs]) ) }))
      }
      L <- BL[,(d+1):(d+K)]
      if (d == 0) {
        Mu.B <- rep(0, p)
      } else {
        Mu.B <- BL[,1:d] %*% t(Cov)
      }
      
      ##Update C##
      C <- t(sapply(seq_len(n), function(i){y <- Y[,i]; ind.obs <- !is.na(y); y <- y[ind.obs] - Mu.B[ind.obs,i]; Li <- L[ind.obs,]; Sigma.i <- Sigma[ind.obs]; solve(t(Li/Sigma.i)%*%Li,t(Li/Sigma.i)%*%y)}))
      W <- cbind(Cov,C)
      
      ##Compute log-likelihood##
      Mu <- BL %*% t(W)
      log.like <- -sum(unlist(lapply(seq_len(p), function(g){1/Sigma[g]*sum((Y[g,!is.na(Y[g,])]-Mu[g,!is.na(Y[g,])])^2)})))
      log.like.vec[j] <- log.like
      
      if (abs(log.like - log.like.0)/abs(log.like.0) < tol) {
        break
      }
      log.like.0 <- log.like
    }
  }
  C <- Q %*% t(Q) %*% C; C <- C %*% solve(chol.default(1/n*t(C)%*%C))
  W <- cbind(Cov,C)
  BL[Frac.missing==0,] <- Y[Frac.missing==0,]%*%W%*%solve(t(W)%*%W)
  
  if (sum(Frac.missing>0) > 0) {
    BL[Frac.missing > 0,] <- t(apply(Y[Frac.missing > 0,], 1, function(y) {ind.obs <- !is.na(y); as.vector( solve(t(W[ind.obs,])%*%W[ind.obs,], t(W[ind.obs,])%*%y[ind.obs]) ) }))
  }
  Sigma <- apply(X = Y - BL%*%t(W), MARGIN = 1, function(x){tmp <- !is.na(x); sum(x[tmp]^2)/(sum(tmp)-d-K)})
  C <- C %*% svd(t(BL[,(d+1):(d+K)]/Sigma)%*%BL[,(d+1):(d+K)])$u
  return(list(C=C, W=cbind(Cov,C), LtSigmaL=svd(t(BL[,(d+1):(d+K)]/Sigma)%*%BL[,(d+1):(d+K)])$d/p, log.like=log.like.vec, Sigma=Sigma))
}

Compute.Q <- function(X) {
  X <- cbind(X)
  qr.X <- qr(X)
  return( qr.Q(qr.X, complete = T)[,(qr.X$rank+1):nrow(X)] )
}
