#####Simulate metabolite data matrix with confounding and missing data, such that residuals are correlated in networks.#####
library(sn)

#Assumes data are (skew) normally distributed
SimMetab.Matrix.conf.corr <- function(n, p, X, Z=NULL, Omega=NULL, log.a.mean=log(1), log.a.sd=0.3, y0.mean=16, y0.sd=1, mu.mean=19, mu.sd=1, Info.vec, sd.L=0.5, mean.rv=1, sd.rv=0.5, 
                                 skew.resid=0, miss.mech=c("logit", "probit", "skew-norm", "t"), skew.miss=0, t.df.miss=7, delta.B=0.05, sigma.B=0.4, mean.B=0, confound.exp=1,
                                 cor.params=list(rho=0,size.block=3,arbitrary.corr=T,n.blocks=NULL,s.rho=0.6)) {
  K <- length(Info.vec)
  out <- list()
  out$X <- X
  X <- cbind(X)
  
  if (is.null(Z)) {Z <- cbind(rep(1,n))}
  out$Z <- Z
  X.tilde <- cbind(Compute.Orthog(X=Z, Z=X))
  
  #Simulate missingness mechanism paramters
  miss.mech <- match.arg(miss.mech, c("logit", "probit", "skew-norm", "t"))
  if (miss.mech == "logit") {log.a.mean <- log.a.mean + log(pi) - 1/2*log(3)}
  if (miss.mech == "t") {log.a.mean <- log.a.mean + 1/2*log(t.df.miss/(t.df.miss-2))}
  out$a <- exp( log.a.mean + log.a.sd*rnorm(p) )
  out$y0 <- y0.mean + y0.sd*rnorm(p)
  if (miss.mech == "skew-norm") {
    if (length(skew.miss) == 1) {skew.miss <- rep(skew.miss, p)}
    out$s <- skew.miss
    out$a <- out$a * sqrt(1-2*skew.miss^2/pi/(1+skew.miss^2))
    out$y0 <- out$y0 - 1/out$a*unlist(lapply(X = skew.miss, function(s){qsn(p = 1/2, alpha = s)}))
  }
  
  #Simulate residual variance#
  if (sd.rv == 0) {
    out$Sigma <- rep(mean.rv, p)
  } else {
    out$Sigma <- rgamma(n = p, shape = mean.rv^2/sd.rv^2, rate = mean.rv/sd.rv^2)
  }
  
  #Simulate mu, L an C
  out$mu <- mu.mean + mu.sd*rnorm(p)
  out$L <- matrix(0, nrow=p, ncol=K)
  out$C <- Simulate.C.conf(n=n, K=K, X=X, Z=Z, Omega=Omega, confound.exp=confound.exp)
  out$Omega <- solve(t(X.tilde)%*%X.tilde, t(X.tilde)%*%out$C)
  out$Var.X <- 1 - sum(X * Compute.Orthog(X=cbind(out$C,Z), Z=X)) / sum(X.tilde^2)
  tt <- matrix(0, nrow=p, ncol=K)
  out$delta.L <- rep(NA, K)
  for (k in 1:K) {
    sd.k <- max(sd.L, sqrt(Info.vec[k]))
    out$delta.L[k] <- min(Info.vec[k]/sd.k^2, 1)
    tmp.ind <- rbinom(n = p, size = 1, prob = out$delta.L[k])
    out$L[tmp.ind==1,k] <- sd.k*rnorm(sum(tmp.ind))*sqrt(out$Sigma[tmp.ind==1])
  }
  out$L <- out$L %*% svd(t(out$L/out$Sigma)%*%out$L)$u
  
  #Simulate B#
  ind.nonzero <- rbinom(n = p, size = 1, prob = delta.B)
  out$B <- ind.nonzero
  ind.neg <- rbinom(n=sum(ind.nonzero==1), size=1, prob=1/2)
  out$B[ind.nonzero == 1][ind.neg == 1] <- (sigma.B * rnorm(sum(ind.neg == 1)) - mean.B) * sqrt(out$Sigma[ind.nonzero == 1][ind.neg == 1])
  out$B[ind.nonzero == 1][ind.neg == 0] <- (sigma.B * rnorm(sum(ind.neg == 0)) + mean.B) * sqrt(out$Sigma[ind.nonzero == 1][ind.neg == 0])
  
  #Simulate Y#
  if (length(skew.resid) == 1 && skew.resid[1] == 0) {
    tmp <- Add.Corr(p = p, n = n, rho = cor.params$rho, size.block = cor.params$size.block, arbitrary.corr = cor.params$arbitrary.corr, t.df = NULL, n.blocks = cor.params$n.blocks, s = cor.params$s.rho)
    out$Y <- out$mu%*%rbind(rep(1,n)) + cbind(out$B)%*%t(X) + out$L%*%t(out$C) + tmp$E*sqrt(out$Sigma)
    out$kappa2 <- tmp$kappa2
  } else {
    if (length(skew.resid) == 1) {
      delta.skew <- skew.resid/sqrt(1+skew.resid^2)
      out$Y <- (out$mu-delta.skew*sqrt(2/pi))%*%rbind(rep(1,n)) + cbind(out$B)%*%t(X) + out$L%*%t(out$C) + 1/sqrt(1-2*delta.skew^2/pi)*matrix(rsn(n = n*p, alpha = skew.resid), nrow=p, ncol=n)*sqrt(out$Sigma)
    } else {
      Resids <- sapply(X = 1:p, function(g){
        s <- skew.resid[g]; delta.s <- s/sqrt(1+s^2)
        return( out$mu[g]-delta.s*sqrt(2/pi) + 1/sqrt(1-2*delta.s^2/pi)*sqrt(out$Sigma[g])*rsn(n = n, alpha = s) )
      })
      out$Y <-  out$mu%*%rbind(rep(1,n)) + cbind(out$B)%*%t(X) + out$L%*%t(out$C) + t(Resids)
    }
  }
  
  #Simulate missing data#
  out$Y.all <- out$Y
  if (miss.mech == "logit") {
    out$Prob.Obs <- t(sapply( X = 1:p, function(g){1/(1+exp( -out$a[g]*(out$Y[g,]-out$y0[g]) ))} ))
  }
  if (miss.mech == "probit") {
    out$Prob.Obs <- t(sapply( X = 1:p, function(g){pnorm(out$a[g]*(out$Y[g,]-out$y0[g]))} ))
  }
  if (miss.mech == "skew-norm") {
    if (length(skew.miss)==1) {skew.miss <- rep(skew.miss[1],p)}
    out$Prob.Obs <- t(sapply( X = 1:p, function(g){psn(x = out$a[g]*(out$Y[g,]-out$y0[g]), alpha = skew.miss[g])} ))
  } 
  if (miss.mech == "t") {
    out$Prob.Obs <- t(sapply( X = 1:p, function(g){pt(out$a[g]*(out$Y[g,]-out$y0[g]), df=t.df.miss)} ))
  }
  out$Y <- t(sapply(X = 1:p, function(g){
    out.g <- out$Y[g,]; out.g[rbinom(n = n, size = 1, prob = out$Prob.Obs[g,]) == 0] <- NA
    return(out.g)
  }))
  out$Frac.Obs <- apply(X = out$Y, MARGIN = 1, function(x){mean(!is.na(x))})
  return(out)
}


####Simulated network residuals####
Add.Corr <- function(p, n, rho, size.block, arbitrary.corr, t.df=NULL, n.blocks, s) {
  out <- list()
  if (is.null(n.blocks)) {n.blocks <- p/size.block}
  if (is.null(t.df)) {
    out$E <- matrix(rnorm(n*p), nrow = p, ncol = n)
  } else {
    out$E <- matrix(rt(n = n*p, df = t.df)*sqrt((t.df-2)/t.df), nrow = p, ncol = n)
  }
  if (!arbitrary.corr) {
    tmp <- matrix(rho, nrow=size.block, ncol=size.block); diag(tmp) <- rep(1, length=size.block)
    R <- CorrConf::sqrt.mat(tmp)
    for (i in 1:n.blocks) {
      out$E[((i-1)*size.block+1):(i*size.block),] <- R%*%out$E[((i-1)*size.block+1):(i*size.block),]
    }
  } else {
    out$kappa2 <- rep(NA, n.blocks)
    for (i in 1:n.blocks) {
      #tmp <- matrix(runif(n = 1, min = min.rho, max = max.rho), nrow=size.block, ncol=size.block); diag(tmp) <- rep(1, length=size.block)
      tmp <- Compute.Corr(s = s, size.block = size.block)
      out$kappa2[i] <- Kappa2(tmp)
      out$E[((i-1)*size.block+1):(i*size.block),] <- (CorrConf::sqrt.mat(tmp))%*%out$E[((i-1)*size.block+1):(i*size.block),]
    }
  }
  return(out)
}

Compute.Corr <- function(s, size.block) {
  go <- 1
  while (go) {
    zz <- s*rnorm(size.block^2); tmp <- matrix((exp(zz)-1)/(exp(zz)+1), nrow=size.block, ncol=size.block); tmp <- 1/2*(tmp + t(tmp)); diag(tmp) <- rep(1, length=size.block)
    go <- as.numeric( min(eigen(tmp)$values) <= 0 )
  }
  return(tmp)
}

Kappa2 <- function(M) {
  t <- eigen(M, symmetric = T)
  return(t$values[1]/t$values[length(t$values)])
}
