#library(quantreg)
#library(sn)
#library(pbivnorm)

###Get starting points for Gamma and C###
#Y is p x n
#Cov is n x d
Theta.Start <- function(Y, Cov, K, regression="ols", parallel=F) {
  iter.svd <- 3
  n <- ncol(Y)
  p <- nrow(Y)
  d <- ncol(Cov)
  p.missing <- apply(Y, 1, function(x){sum(is.na(x))})
  Sigma <- rep(NA, p)
  Gamma <- matrix(NA, nrow=p, ncol=d+K)
  
  #Get starting point for C#
  if (K > 0) {
    Q <- qr.Q(qr(Cov), complete = T)[,(d+1):n]
    Y.0 <- Y[p.missing == 0,] %*% Q
    for (i in 1:iter.svd) {
      if (i == 1) {
        C <- sqrt(n) * cbind(svd(Y.0)$v[,1:K])
        Sigma.0 <- Est.Sigma.0(Y=Y.0, Cov=C)
      } else {
        C <- sqrt(n) * cbind(svd(Y.0 / sqrt(Sigma.0))$v[,1:K])
        Sigma.0 <- Est.Sigma.0(Y=Y.0, Cov=C)
      }
    }
    C <- Q%*%C
    Sigma[p.missing == 0] <- Sigma.0
    Cov.total <- cbind(Cov,C)
  } else {
    Cov.total <- cbind(Cov)
    C <- NULL
    Sigma[p.missing == 0] <- Est.Sigma.0(Y=Y[p.missing == 0,], Cov=Cov.total)
  }
  Gamma[p.missing == 0,] <- Y[p.missing == 0,] %*% Cov.total %*% solve(t(Cov.total)%*%Cov.total)
  if (regression != "ols") {
    Y.tmp <- Y[p.missing > 0,]; Y.tmp[is.na(Y.tmp)] <- 0
    Gamma[p.missing > 0,] <- Med.Regression(Y=Y.tmp, Cov=Cov.total, parallel=parallel)
    rm(Y.tmp)
  } else {
    tmp <- t(apply(rbind(Y[p.missing > 0,]), 1, function(y){ ind.use <- !is.na(y); solve(t(Cov.total[ind.use,])%*%Cov.total[ind.use,], t(Cov.total[ind.use,])%*%y[ind.use]) }))
    if (ncol(tmp) != ncol(Gamma)) {tmp <- t(tmp)}
    Gamma[p.missing > 0,] <- tmp
  }
  Sigma[p.missing > 0] <- Est.Sigma.0(Y = rbind(Y[p.missing > 0,]), Cov=Cov.total, Gamma=Gamma[p.missing > 0,])
  
  return(list( C=C, W=Cov.total, Gamma=Gamma, Sigma=Sigma ))
}


############   Estimate y0 and a   ############

Est.ay0 <- function(a, y.0, Y, Mu, Sigma, Pi, min.y0 = 1, min.a0 = 0.01, max.iter=100, tol=1e-6) {
  out <- constrOptim(theta=c(a,y.0), mlike.theta, mgrad.theta, ui = diag(2), ci = c(min.a0, min.y0), method = "BFGS", control = list(maxit=max.iter, reltol=tol), Y=Y, Mu=Mu, Sigma=Sigma, Pi=Pi)
  return(list(a=out$par[1], y.0=out$par[2], out=out$convergence))
}

mlike.theta <- function(theta, Y, Mu, Sigma, Pi) {
  a <- theta[1]
  y.0 <- theta[2]
  p <- nrow(Y)
  n <- ncol(Y)

  Ind <- matrix(1, nrow=p, ncol=n); Ind[is.na(Y)] <- 0 
  ind.obs <- matrix(Ind == 1, nrow=p, ncol=n)
  ind.unobs <- matrix(Ind == 0, nrow=p, ncol=n)
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  Pi <- cbind(Pi) %*% rbind(rep(1,n))
  
  Quant <- a * Y - y.0; Quant[ind.unobs] <- -(a*Mu[ind.unobs] - y.0) / sqrt(1+a^2*Sigma[ind.unobs])
  like <- sum( Pi[ind.obs] * pnorm(Quant[ind.obs], log.p = T) ) + sum( Pi[ind.unobs] * pnorm(Quant[ind.unobs], log.p = T) )
  return(-like)
}

mgrad.theta <- function(theta, Y, Mu, Sigma, Pi) {
  a <- theta[1]
  y.0 <- theta[2]
  p <- nrow(Y)
  n <- ncol(Y)
  
  Ind <- matrix(1, nrow=p, ncol=n); Ind[is.na(Y)] <- 0 
  ind.obs <- matrix(Ind == 1, nrow=p, ncol=n)
  ind.unobs <- matrix(Ind == 0, nrow=p, ncol=n)
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  Pi <- cbind(Pi) %*% rbind(rep(1,n))
  
  Cmat <- a * Y - y.0; Cmat[ind.unobs] <- (a * Mu[ind.unobs] - y.0) / sqrt(1 + a^2*Sigma[ind.unobs])
  phi.d.Phi <- matrix(0,nrow=p,ncol=n); phi.d.Phi[ind.obs] <- dnorm(Cmat[ind.obs])/pnorm(Cmat[ind.obs]); phi.d.Phi[ind.unobs] <- dnorm(Cmat[ind.unobs])/pnorm(-Cmat[ind.unobs])
  
  grad.a <- sum( Pi[ind.obs] * Y[ind.obs] * phi.d.Phi[ind.obs] ) - sum( Pi[ind.unobs]*( Mu[ind.unobs] + a*Sigma[ind.unobs]*y.0 )/(1+a^2*Sigma[ind.unobs])^1.5 * phi.d.Phi[ind.unobs] )
  grad.y0 <- -sum(Pi[ind.obs] * phi.d.Phi[ind.obs]) + sum(Pi[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs])*phi.d.Phi[ind.unobs])
  return(-c(grad.a,grad.y0))
}


#######   Estimate posterior probs   #######
Est.Post <- function(q, a, Y, Mu, y.0, Sigma, prior.nonrand) {
  p <- nrow(Y)
  n <- ncol(Y)
  if (prior.nonrand > (1-1e-8)) {return(rep(1,p))}
  if (prior.nonrand < 1e-8) {return(rep(0,p))}
  
  Ind <- matrix(1, nrow=p, ncol=n); Ind[is.na(Y)] <- 0
  
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  ind.obs <- Ind == 1
  ind.unobs <- Ind == 0
  N.vec <- apply(Ind, 1, sum)
  
  c.mat <- a*Y - y.0
  c.mat[ind.unobs] <- (a*Mu[ind.unobs] - y.0) / sqrt(1 + a^2*Sigma[ind.unobs])
  log.Probs <- matrix(NA, nrow=p, ncol=n)
  log.Probs[ind.obs] <- pnorm(c.mat[ind.obs], log.p = T); log.Probs[ind.unobs] <- pnorm(-c.mat[ind.unobs], log.p = T)
  
  unlist( lapply(seq_len(p), function(g){ l.nonr <- sum(log.Probs[g,]); l.r <- (n-N.vec[g])*log(1-q)+N.vec[g]*log(q); C.max <- max(l.nonr,l.r); l.nonr <- l.nonr - C.max; l.r <- l.r - C.max; return(prior.nonrand*exp(l.nonr)/(prior.nonrand*exp(l.nonr) + (1-prior.nonrand)*exp(l.r))) }) )
}

######    Estimate q and pi    ######
Est.q <- function(Pi, Y) {
  N.vec <- apply(Y, 1, function(x){sum(!is.na(x))})
  return( sum( (1-Pi)*N.vec/ncol(Y)/sum(1-Pi) ) )
}

Est.pi <- function(Pi) {
  mean(Pi)
}


#########   Update Mu   #########
Est.Mu <- function(C, Gamma, Cov, Y, Sigma, a, y.0, Pi, max.iter=100) {
  p <- nrow(Y)
  n <- ncol(Y)
  Ind <- matrix(0,nrow=p,ncol=n)
  Ind[!is.na(Y)] <- 1

  for (i in 1:max.iter) {
    Cov.total <- cbind(Cov, C)
    Mu <- cbind(Gamma) %*% rbind(t(Cov.total))
    Resids <- Y - Mu
    for (g in 1:p) {
      ind.g <- Ind[g,]
      v.g <- Sigma[g]
      C.g <- (a*Mu[ind.g==0] - y.0)/sqrt(1+a^2*v.g)
      
      x <- rep(NA,g)
      x[ind.g==1] <- Resids[g,ind.g==1]/v.g
      x[ind.g==0] <- -Pi[g]*a^2/sqrt(1+a^2*v.g)*(dnorm(C.g)/pnorm(-C.g))
      grad.g <- t(Cov.total) %*% x
    }
  }
}


###Median regression###
#Y is p x n
#Cov is n x d
Med.Regression <- function(Y, Cov, parallel=F) {
  p <- nrow(Y)
  Y.list <- lapply(seq_len(p), function(g){Y[g,]})
  rm(Y)
  
  if (parallel) {
    n_cores <- max(detectCores() - 1, 1)
    cl <- makeCluster(n_cores)
    clusterExport(cl, "Cov")
    clusterEvalQ(cl, library(quantreg))
    out <- t(parSapply(cl=cl, Y.list, function(y, Cov) { suppressWarnings(rq.fit(x=Cov, y=y)$coefficients) }, Cov=Cov))
    stopCluster(cl)
  } else {
    out <- t(sapply(Y.list, function(y) {suppressWarnings(rq.fit(x=Cov, y=y)$coefficients)}))
  }
  
  return(out)
}


###Estimate Sigma###
Est.Sigma.0 <- function(Y, Cov, Gamma=NULL) {
  Sigma <- rep(NA, nrow(Y))
  p.missing <- apply(Y,1,function(x){sum(is.na(x))})
  if (sum(p.missing == 0) > 0) {
    R <- rbind(Y[p.missing == 0,] - Y[p.missing == 0,] %*% Cov %*% solve(t(Cov)%*%Cov,t(Cov)))
    Sigma[p.missing == 0] <- rowSums(R^2)/(nrow(Cov)-ncol(Cov))
  }
  if (sum(p.missing > 0) == 0) {return(Sigma)}
  Gamma <- rbind(Gamma)
  if (nrow(Gamma) != nrow(Y)) {Gamma <- t(Gamma)}
  Sigma[p.missing > 0] <- apply(Y[p.missing > 0,] - Gamma[p.missing > 0,] %*% t(Cov), 1, function(y){ sum(y[!is.na(y)]^2)/sum(!is.na(y)) })
  return(Sigma)
}

dnorm.prime <- function(x) {
  return(-x/sqrt(2*pi) * exp(-x^2/2))
}



###################   Skew normal missingness   ###################

Est.skew.miss <- function(theta=c(1,15,1), Y, Mu, Sigma, lambda=1, K.ind=NULL, C=NULL, max.ignore = 0.05, s=NULL, max.s=5, min.s=-5, min.a=1e-8, min.y0=-5) {
  Prob.missing <- apply(X = Y, MARGIN = 1, function(x){sum(is.na(x))})/ncol(Y)
  ind.use <- Prob.missing == 0 | Prob.missing > max.ignore
  if (!is.null(s)){
    out <- constrOptim(theta = theta[1:2], f = mlike.miss, grad = mgrad.miss, ui = diag(2), ci = c(min.a, min.y0), Y=rbind(Y[ind.use,]), Mu=rbind(Mu[ind.use,]), Sigma=Sigma[ind.use], s=s)
    return(list(a=out$par[1], y.0=out$par[2], s=s, out=out$convergence, mlike=out$value))
  }
  
  if (is.null(C)) {
    U <- cbind(rep(1,ncol(Y)))
    W = NULL
  } else {
    U <- Get.U(C = C, K.ind = K.ind)
    W <- 1/nrow(U)*t(U)%*%U
  }
  out <- constrOptim(theta = theta, f = mlike.miss, grad = mgrad.miss, ui = rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(0,0,-1)), ci = c(min.a, min.y0, min.s, -max.s), Y=rbind(Y[ind.use,]), Mu=rbind(Mu[ind.use,]), Sigma=Sigma[ind.use], lambda=lambda, U=U, W=W)
  return(list(a=out$par[1], y.0=out$par[2], s=out$par[3], out=out$convergence, mlike=out$value))
}

#Minus log-likelihood for missingness
mlike.miss <- function(theta, Y, Mu, Sigma, s=NULL, lambda, U, W=NULL) {
  a <- theta[1]; y.0 <- theta[2]
  if (is.null(s)) {
    s <- theta[3]
  }
  ind.obs <- !is.na(Y)
  prob.mat <- t(sapply(seq_len(nrow(Y)), function(g, a, y.0, s) { Compute.Prob.skew(y=Y[g,], mu=Mu[g,], v=Sigma[g], a=a, y.0=y.0, s=s) }, a=a, y.0=y.0, s=s))
  prob.mat[prob.mat <= 0] <- 1e-8
  mlog.like <- -(sum(log(prob.mat[ind.obs])) + sum(log(1-prob.mat[!ind.obs])))
  
  if (lambda < 1) {
    g.GMM <- 2*func.GMM(theta = theta, y = Y[1,], U = U, W = W)
  } else {
    g.GMM <- 0
  }
  return(lambda*mlog.like + (1-lambda)*t(g.GMM)%*%g.GMM)
}

#Minus gradient of log-likelihood for missingness
mgrad.miss <- function(theta, Y, Mu, Sigma, s=NULL, lambda, U, W=NULL) {
  a <- theta[1]; y.0 <- theta[2]
  if (is.null(s)) {
    choose.s <- 1
    s <- theta[3]
  } else {
    choose.s <- 0
  }
  n <- ncol(Y)
  p <- nrow(Y)
  ind.obs <- !is.na(Y)
  ind.unobs <- !ind.obs
  
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  C <- a*Y - y.0; C[ind.unobs] <- (a*Mu[ind.unobs] - y.0)/sqrt(1+a^2*Sigma[ind.unobs])
  
  prob.mat <- t(sapply(seq_len(nrow(Y)), function(g, a, y.0, s) { Compute.Prob.skew(y=Y[g,], mu=Mu[g,], v=Sigma[g], a=a, y.0=y.0, s=s) }, a=a, y.0=y.0, s=s))
  prob.mat[prob.mat <= 0] <- 1e-8
  
  #skew normal for observed data#
  d.skew <- dsn(C[ind.obs], alpha = s)
  
  #Initialize gradient with observed data#
  dl.da <- sum(Y[ind.obs] * d.skew/prob.mat[ind.obs])
  dl.dy0 <- -sum(d.skew/prob.mat[ind.obs])
  dl.ds <- -sqrt(2/pi)/(1+s^2) * sum(dnorm(C[ind.obs]*sqrt(1+s^2))/prob.mat[ind.obs])
  
  #Gradient of f w.r.t a#
  df.da.1 <- (y.0 + C[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs])) * pnorm(s*C[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]))
  df.da.2 <- s*a^2*Sigma[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs])/sqrt(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]) * dnorm(s*C[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]))
  df.da <- dnorm(C[ind.unobs])/a/sqrt(1+a^2*Sigma[ind.unobs]) * ( df.da.1 + df.da.2 )
  
  #Gradient of f w.r.t y.0#
  df.dy0 <- -dnorm(C[ind.unobs])/sqrt(1+a^2*Sigma[ind.unobs]) * pnorm(s*C[ind.unobs]/sqrt(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]))
  
  #Gradient of f w.r.t s#
  if (choose.s) {
    df.ds <- -1/sqrt(2*pi)/(s^2+1)/sqrt(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]) * dnorm(C[ind.unobs]*sqrt( (s^2+1)*(1+a^2*Sigma[ind.unobs])/(1+a^2*Sigma[ind.unobs]+s^2*a^2*Sigma[ind.unobs]) ))
    mgrad.loglike <- -c( dl.da - 2*sum(df.da/(1-prob.mat[ind.unobs])), dl.dy0 - 2*sum(df.dy0/(1-prob.mat[ind.unobs])), dl.ds - 2*sum(df.ds/(1-prob.mat[ind.unobs])) )
    if (lambda < 1) {
      grad.GMM <- 2*Grad.GMM(theta = theta, y = Y[1,], U = U, W = W)
    } else {
      grad.GMM <- 0
    }
    return(lambda*mgrad.loglike + (1-lambda)*grad.GMM)
    return(  )
  } else {
    return( -c( dl.da - 2*sum(df.da/(1-prob.mat[ind.unobs])), dl.dy0 - 2*sum(df.dy0/(1-prob.mat[ind.unobs])) ) )
  }
}

Compute.Prob.skew <- function(y, mu, v, a, y.0, s) {  #Probability metabolite is observed
  n <- length(y)
  out <- rep(NA, n)
  obs.y <- !is.na(y)
  out[obs.y] <- psn(x=a*y[obs.y]-y.0, alpha = s)
  if (sum(obs.y) == n) {return(out)}
  d <- 1/a/sqrt(v)
  out[!obs.y] <- 2*pbivnorm(x = 0, y = d*(a*mu[!obs.y] - y.0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE)
  return(out)
}

Prob.skew <- function(a, y.0, s, Y, Mu, Sigma) {
  return( t( sapply(seq_len(nrow(Y)), function(g){Compute.Prob.skew(y=Y[g,], mu=Mu[g,], v=Sigma[g], a=a, y.0=y.0, s=s)}) ) )
}


###########   Inverse probability weighting to estimate L, B, C and Sigma   ###########

IPW.Est <- function(Y, Cov, Pi, C, max.ignore=0.05, max.iter=100, like.tol=1e-5) {
  p <- nrow(Y)
  n <- ncol(Y)
  d <- ncol(Cov)
  K <- ncol(C)
  
  Ind.obs <- matrix(1, nrow=p, ncol=n); Ind.obs[is.na(Y)] <- 0
  Prop.missing <- 1 - rowSums(Ind.obs)/n
  Missing <- which(Prop.missing > 0)
 #Pi[Prop.missing<=max.ignore,] <- 1
  N.eff <- rowSums(1/Pi * Ind.obs)
  
  C.total <- cbind(Cov,C)
  like <- -1e16
  
  for (i in 1:max.iter) {
    #Update [B L]#
    BL <- matrix(NA, nrow=p, ncol=d+K)
    BL[-Missing,] <- t(solve(t(C.total)%*%C.total,t(C.total)%*%t(Y[-Missing,])))
    BL[Missing,] <- t( sapply(X = Missing, function(g){ ind.g <- Ind.obs[g,]==1; pi.g <- Pi[g,ind.g]; solve( t(C.total[ind.g,] / pi.g)%*%C.total[ind.g,], t(C.total[ind.g,] / pi.g) %*% Y[g,ind.g] ) }) )
    B <- cbind(BL[,1:d])
    L <- cbind(BL[,(d+1):(d+K)])
    
    #Update C#
    C <- t( sapply(X = seq_len(n), function(i){ ind.i <- Ind.obs[,i]==1; pi.i <- Pi[ind.i,i]; solve( t(L[ind.i,] / pi.i)%*%L[ind.i,], t(L[ind.i,] / pi.i)%*%(Y[ind.i,i] - B[ind.i,]%*%cbind(Cov[i,])) ) }) )
    C.total <- cbind(Cov,C)
    Mu <- BL %*% t(C.total)
    
    #Update Sigma#
    Sigma <- rowSums((Y - Mu)^2 / Pi, na.rm = T) / N.eff
    
    #Compute weighted likelihood#
    like.i <- -sum(N.eff * log(Sigma)) - sum( ((Y - Mu)^2/Sigma)/Pi, na.rm = T )
    if (abs(like.i - like)/abs(like) <= like.tol) {
      return(list( C.total=C.total, BL=BL, Sigma=Sigma, n.iter=i, out=0 ))
    }
    like <- like.i
  }
  return(list( C.total=C.total, BL=BL, Sigma=Sigma, n.iter=i, out=1 ))
}


###################   Junk   ###################

########    Estmiate y_0    ########
Est.y0 <- function(y0, Y, Mu, Sigma, a, Pi, max.iter=100, tol=1e-5) {
  p <- nrow(Y)
  n <- ncol(Y)
  p.nonrandom <- sum(Pi)
  
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  Pi <- cbind(Pi) %*% rbind(rep(1,n))
  Ind <- matrix(1, nrow=p, ncol=n); Ind[is.na(Y)] <- 0
  ind.obs <- Ind == 1
  ind.unobs <- Ind == 0
  
  for (i in 1:max.iter) {
    Delta <- Y - y0; Delta[ind.unobs] <- Mu[ind.unobs] - y0
    
    #Create c#
    c.mat <- a * Delta
    c.mat[ind.unobs] <- c.mat[ind.unobs] / sqrt(1 + a^2*Sigma[ind.unobs])
    dens.c <- dnorm(c.mat)
    dens.prime.c <- dnorm.prime(c.mat)
    prob.c <- pnorm(c.mat)
    
    #gradient#
    grad <- -a * sum(Pi[ind.obs] * dens.c[ind.obs]/prob.c[ind.obs]) + a*sum(Pi[ind.unobs] / sqrt(1+a^2*Sigma[ind.unobs]) * dens.c[ind.unobs]/(1-prob.c[ind.unobs]))
    if (abs(grad)/p.nonrandom < tol) {
      return(list(y0=y0, out=0, n.iter=i))
    }
    
    #hessian#
    hess <- a^2 * sum(Pi[ind.obs] * ( dens.prime.c[ind.obs]/prob.c[ind.obs] - dens.c[ind.obs]^2/prob.c[ind.obs]^2 ))
    hess <- hess - a^2 * sum(Pi[ind.unobs] / (1+a^2*Sigma[ind.unobs]) * ( dens.prime.c[ind.unobs]/(1-prob.c[ind.unobs]) + dens.c[ind.unobs]^2/(1-prob.c[ind.unobs])^2 ))
    
    y0 <- y0 - grad/hess
  }
  return(list(y0=y0, out=1, n.iter=i))
}


############   Estimate a   ############
#Delta_{gi} = y_{gi} - y_0 if observed (i.e. Ind_{gi} = 1) and mu_{gi} - y_0 if Ind_{gi} = 0
#Pi as a p-vector of posterior probability for the missingness mechanism for the metabolite
#Sigma is a p-vector of residual variances
Est.a <- function(a.0, Y, Mu, y0, Sigma, Pi, max.iter=100, tol=1e-5, c.wolfe=1/2) {
  p <- nrow(Y)
  n <- ncol(Y)
  p.nonrandom <- sum(Pi)
  
  Delta <- Y - y0; Delta[is.na(Y)] <- Mu[is.na(Y)] - y0
  Ind <- matrix(1, nrow=p, ncol=n); Ind[is.na(Y)] <- 0
  
  Sigma <- cbind(Sigma) %*% rbind(rep(1,n))
  Pi <- cbind(Pi) %*% rbind(rep(1,n))
  ind.obs <- matrix(Ind == 1, nrow=p, ncol=n)
  ind.unobs <- matrix(Ind == 0, nrow=p, ncol=n)
  
  for (i in 1:max.iter) {
    #Create c#
    c.mat <- a.0 * Delta
    c.mat[ind.unobs] <- c.mat[ind.unobs] / sqrt(1 + a.0^2*Sigma[ind.unobs])
    dens.c <- dnorm(c.mat)
    prob.c <- pnorm(c.mat)
    
    #Gradient#
    grad <- sum( Pi[ind.obs] * Delta[ind.obs] * dens.c[ind.obs] / prob.c[ind.obs] ) - sum( Pi[ind.unobs] * Delta[ind.unobs] / (1 + a.0^2*Sigma[ind.unobs])^1.5 * dens.c[ind.unobs] / (1 - prob.c[ind.unobs]) )
    if (abs(grad)/p.nonrandom < tol) {return(list(a=a.0, out=0, n.iter=i))}
    
    #Update a#
    like.0 <- like.a(a.0, Delta, Sigma, Pi, Ind)
    for (j in 1:100) {
      step.j <- sign(grad) * c.wolfe^(j-1)
      a.test <- max(step.j + a.0, 0.001)
      if (like.a(a.test, Delta, Sigma, Pi, Ind) > like.0) {
        a <- a.test
        break
      }
    }
    if (abs(a - a.0) < 0.001) {return(list(a=a, out=0.5, n.iter=i))}
    a.0 <- a
  }
  return(list(a=a.0, out=1, n.iter=i))
}

like.a <- function(a, Delta, Sigma, Pi, Ind) {
  p <- nrow(Delta)
  n <- ncol(Delta)
  
  ind.obs <- Ind == 1
  ind.unobs <- Ind == 0
  
  c.mat <- a * Delta
  c.mat[ind.unobs] <- c.mat[ind.unobs] / sqrt(1 + a^2*Sigma[ind.unobs])
  
  return( sum(Pi[ind.unobs] * pnorm(-c.mat[ind.unobs], log.p = T)) + sum(Pi[ind.obs] * pnorm(c.mat[ind.obs], log.p = T)) )
}



############ Estimate L, C and Sigma from the full data likelihood, given missingness parameters ############
#Gamma = [Beta, L]
EstParam.Full <- function(Y, Cov=NULL, C=NULL, Gamma, Sigma, theta.missing, min.missing.perc=0.05, max.iter=100, like.tol=1e-5) {
  Gamma <- cbind(Gamma)
  Cov <- cbind(Cov)
  C <- cbind(C)
  
  p <- nrow(Y)
  n <- ncol(Y)
  K <- ncol(C)
  d <- ncol(Cov)
  
  Prob.miss.true <- apply(Y,1,function(x){sum(is.na(x))})/ncol(Y)
  
  if (!is.null(C)) {L <- Gamma[,(d+1):(K+d)]}
  if (!is.null(Cov)) {Beta <- Gamma[,1:d]}
  Mu <- Gamma %*% t(cbind(Cov,C))
  log.like <- -1e16*p*n
  for (i in 1:max.iter) {
    ##Update Gamma##
    Gamma <- cbind(t(sapply(X = seq_len(p), function(g){ if (is.na(theta.missing[g,1])) { tmp.y <- Y[g,]; ind.obs <- !is.na(tmp.y); tmp.cov <- cbind(Cov,C)[ind.obs,]; return( solve(t(tmp.cov)%*%tmp.cov, t(tmp.cov)%*%tmp.y[ind.obs]) ) }; out <- optim(par = Gamma[g,], fn = mgamma.like, gr = mgamma.grad, y = Y[g,], Cov = cbind(Cov,C), theta.missing = theta.missing[g,], v=Sigma[g], method = "BFGS"); out$par })))
    
    ##Update C##
    if (!is.null(C)) {
      C <- cbind(t(sapply(X = seq_len(n), function(i){ if (sum(is.na(Y[,i]))==0){ if (is.null(Cov)) {return( solve(t(L/Sigma)%*%L, t(L/Sigma)%*%Y[,i]) )}; return( solve(t(L/Sigma)%*%L, t(L/Sigma)%*%(Y[,i] - Beta%*%cbind(Cov[i,]))) ) }; out <- optim(par = C[i,], fn = mc.like, gr = mc.grad, y = Y[,i], Gamma = Gamma, cov = Cov[i,], theta.missing = theta.missing, Sigma = Sigma, method = "BFGS"); out$par })))
    }
    Mu <- Gamma %*% t(cbind(Cov,C))
    if (!is.null(C)) {L <- cbind(Gamma[,(d+1):(K+d)])}
    if (!is.null(Cov)) {Beta <- cbind(Gamma[,1:d])}
    
    ##Update Sigma##
    Sigma <- as.vector(sapply(X = seq_len(p), function(g){ if (is.na(theta.missing[g,1])) { tmp.y <- Y[g,]; ind.obs <- !is.na(tmp.y); return(sum((tmp.y[ind.obs]-Mu[g,ind.obs])^2)/sum(ind.obs)) }; out <- optimize(f = log.like.v, interval = c(Sigma[g]/10,10*Sigma[g]), maximum = T, y = Y[g,], mu = Mu[g,], theta.missing = theta.missing[g,]); out$maximum }))
  }
  return(list(Gamma=Gamma, C=C, Sigma=Sigma, W=cbind(Cov,C)))
}

##Minus likelihood and gradient for C##

mc.like <- function(par, y, Gamma, cov, theta.missing, Sigma) {
  Gamma <- cbind(Gamma)
  ind.rem <- is.na(y) & is.na(theta.missing[,1])
  if (sum(ind.rem) > 0) {
    theta.missing <- theta.missing[!ind.rem,]
    Gamma <- Gamma[!ind.rem,]
    y <- y[!ind.rem]
    Sigma <- Sigma[!ind.rem]
  }
  
  a <- theta.missing[,1]
  y0 <- theta.missing[,2]
  s <- theta.missing[,3]
  
  mu <- as.vector( Gamma%*%cbind(c(cov,par)) )
  ind.obs <- !is.na(y)
  d <- 1/a[!ind.obs]/sqrt(Sigma[!ind.obs])
  
  log.like <- sum(log( 1 - 2*pbivnorm(x = 0, y = d*(a[!ind.obs]*mu[!ind.obs] - y0[!ind.obs])/sqrt(d^2+1), rho = -s[!ind.obs]*d/sqrt(s[!ind.obs]^2+1)/sqrt(d^2+1), recycle = TRUE) ))
  log.like <- log.like - sum(1/2/Sigma[ind.obs]*(y[ind.obs]-mu[ind.obs])^2)
  return(-log.like)
}

mc.grad <- function(par, y, Gamma, cov, theta.missing, Sigma) {
  Gamma <- cbind(Gamma)
  ind.rem <- is.na(y) & is.na(theta.missing[,1])
  if (sum(ind.rem) > 0) {
    theta.missing <- theta.missing[!ind.rem,]
    Gamma <- Gamma[!ind.rem,]
    y <- y[!ind.rem]
    Sigma <- Sigma[!ind.rem]
  }
  
  a <- theta.missing[,1]
  y0 <- theta.missing[,2]
  s <- theta.missing[,3]
  
  mu <- as.vector( Gamma%*%cbind(c(cov,par)) )
  ind.obs <- !is.na(y)
  d <- 1/a[!ind.obs]/sqrt(Sigma[!ind.obs])
  
  prob.miss <- 1 - 2*pbivnorm(x = 0, y = d*(a[!ind.obs]*mu[!ind.obs] - y0[!ind.obs])/sqrt(d^2+1), rho = -s[!ind.obs]*d/sqrt(s[!ind.obs]^2+1)/sqrt(d^2+1), recycle = TRUE)
  c.vec <- (a[!ind.obs]*mu[!ind.obs]-y0[!ind.obs])/sqrt(1+a[!ind.obs]^2*Sigma[!ind.obs])
  d.vec <- s[!ind.obs]*c.vec/sqrt(1+a[!ind.obs]^2*Sigma[!ind.obs]+s[!ind.obs]^2*a[!ind.obs]^2*Sigma[!ind.obs])
  
  grad.log.like <- as.vector(t(Gamma[ind.obs,]) %*% cbind((y[ind.obs]-mu[ind.obs])/Sigma[ind.obs])) - as.vector(t(Gamma[!ind.obs,]) %*% cbind(1/prob.miss * 2*a[!ind.obs]/sqrt(1+a[!ind.obs]^2*Sigma[!ind.obs]) * dnorm(c.vec) * pnorm(d.vec)))
  return(-grad.log.like[(ncol(Gamma)-length(par)+1):ncol(Gamma)])
}

##Minus likelihood and gradient for Gamma##

mgamma.like <- function(par, y, Cov, theta.missing, v) {
  mu <- as.vector( cbind(Cov)%*%cbind(par) )
  ind.obs <- !is.na(y)
  a <- theta.missing[1]
  y0 <- theta.missing[2]
  s <- theta.missing[3]
  d <- 1/a/sqrt(v)
  
  log.like <- sum(log( 1 - 2*pbivnorm(x = 0, y = d*(a*mu[!ind.obs] - y0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE) ))
  log.like <- log.like - sum(1/2/v*(y[ind.obs]-mu[ind.obs])^2)
  return(-log.like)
}

mgamma.grad <- function(par, y, Cov, theta.missing, v) {
  Cov <- cbind(Cov)
  mu <- as.vector( cbind(Cov)%*%cbind(par) )
  ind.obs <- !is.na(y)
  a <- theta.missing[1]
  y0 <- theta.missing[2]
  s <- theta.missing[3]
  d <- 1/a/sqrt(v)
  
  prob.miss <- 1 - 2*pbivnorm(x = 0, y = d*(a*mu[!ind.obs] - y0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE)
  c.vec <- (a*mu[!ind.obs]-y0)/sqrt(1+a^2*v)
  d.vec <- s*c.vec/sqrt(1+a^2*v+s^2*a^2*v)
  grad.log.like <- as.vector(t(Cov[ind.obs,]) %*% cbind( (y[ind.obs]-mu[ind.obs])/v )) - as.vector(t(Cov[!ind.obs,]) %*% cbind( 1/prob.miss * 2*a/sqrt(1+v*a^2) * dnorm(c.vec) * pnorm(d.vec) ))
  return(-grad.log.like)
}

##Optimize variance##

optim.v <- function(v, y, mu, theta.missing, max.iter=100, tol.grad=1e-6, c.wolfe=0.5) {
  ind.obs <- !is.na(y)
  a <- theta.missing[1]
  y0 <- theta.missing[2]
  s <- theta.missing[3]
  n <- length(y)
  
  resids <- (y[ind.obs]-mu[ind.obs])
  
  for (i in 1:max.iter) {
    grad.i <- 1/2 * sum(resids^2/v^2 - 1/v)
    d <- 1/a/sqrt(v)
    prob.miss <- 1 - 2*pbivnorm(x = 0, y = d*(a*mu[!ind.obs] - y0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE)
    c.vec <- (a*mu[!ind.obs]-y0)/sqrt(1+a^2*v)
    d.vec <- s*c.vec/sqrt(1+a^2*v+s^2*a^2*v)
    grad.i <- grad.i + sum(1/prob.miss * dnorm(c.vec) * ( a^2*c.vec/(1+a^2*v)*pnorm(d.vec) - s*a^2/sqrt(1+a^2*v)/sqrt(1+a^2*v+s^2*a^2*v)*dnorm(d.vec) ))
    if (abs(grad.i/n) <= tol.grad) {
      return(list(par=v, out=0, n.iter=i))
    }
    
    #Update v#
    alpha <- max(sign(grad.i), 1e-8-v)
    count <- 0
    log.like.i <- log.like.v(v, y=y, mu=mu, theta.missing=theta.missing)
    while(log.like.i > log.like.v(v=v+alpha, y=y, mu=mu, theta.missing=theta.missing) & count < 100) {
      alpha <- alpha*c.wolfe
      count <- count + 1
    }
    v <- v + alpha
  }
  
  return(list(par=v, out=1, n.iter=i))
}

log.like.v <- function(v, y, mu, theta.missing) {
  ind.obs <- !is.na(y)
  a <- theta.missing[1]
  y0 <- theta.missing[2]
  s <- theta.missing[3]
  d <- 1/a/sqrt(v)
  
  prob.miss <- 1 - 2*pbivnorm(x = 0, y = d*(a*mu[!ind.obs] - y0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE)
  sum(log(prob.miss)) - sum(1/2*log(v) + 1/2/v*(y[ind.obs]-mu[ind.obs])^2)
}

###Full data log-likelihood###
log.like.full <- function(y, mu, v, theta.missing) {
  ind.obs <- !is.na(y)
  a <- theta.missing[1]
  y0 <- theta.missing[2]
  s <- theta.missing[3]
  d <- 1/a/sqrt(v)
  
  prob.miss <- 1 - 2*pbivnorm(x = 0, y = d*(a*mu[!ind.obs] - y0)/sqrt(d^2+1), rho = -s*d/sqrt(s^2+1)/sqrt(d^2+1), recycle = TRUE)
  prob.obs <- psn(x=a*y[ind.obs]-y0, alpha = s)
  sum(log(prob.miss)) + sum(log(prob.obs)) - sum(1/2*log(v) + 1/2/v*(y[ind.obs]-mu[ind.obs])^2)
}
