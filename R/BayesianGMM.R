####Estimate missingness parameters using a Bayesian model with empirical Bayes priors####
#V.prior is the 2x2 prior variance for (a, y0)
Bayes.GMM <- function(y, C, K.ind, p.min.1=0, p.min.2=0, t.df=4, y0.mean, log.a.mean, V, y0.start=NULL, a.start=NULL, prop.y0.sd=0.2, prop.a.sd=0.2, n.iter=1e4, n.burn=1e3, save.every=10, include.norm=T, min.prob=1/n) {
  U <- Get.U(C = cbind(C), K.ind = K.ind)
  max.store <- 1e4
  U <- cbind(U)
  n <- length(y)
  if (p.min.1 <= 0) {p.min.1 <- 0; p.min.2 <- 0}
  if (p.min.1 > 0 || min.prob <= 0) {min.prob <- 0}
  if (save.every <= 0) {
    save.every <- 0
  } else {
    save.every <- round(save.every)
    if (floor(n.iter/save.every) >= max.store) {save.every <- floor(n.iter/max.store)}
  }
  log.a.sd <- sqrt(V[1,1])
  y0.sd <- sqrt(V[2,2])
  cond.var.log.a <- V[1,1] - V[1,2]^2/V[2,2]
  cond.var.y0 <- V[2,2] - V[1,2]^2/V[1,1]
  ind.obs <- !is.na(y)
  min.y <- min(y[ind.obs])
  
  #Define output#
  out <- list()
  out$Post.Exp <- c(0,0)  #Posterior expectation of (a, y0)
  Post.Exp2 <- matrix(0, nrow=2, ncol=2)  #Posterior expectation of (a, y0)^T(a, y0); a list of 2 x 2 matrices
  out$Post.Exp.W <- rep(0,n)  #Posteior expectations of pi^{-1}
  out$Post.Exp.Pi <- rep(0,n)  #Posteior expectations of pi
  Post.Exp.W2 <- rep(0,n)  #Posteior expectations of pi^{-2}
  
  if (save.every) {
    out$Samples <- list(); n.save <- floor(n.iter/save.every); out$Samples$a <- rep(NA, n.save); out$Samples$y0 <- rep(NA, n.save)
  }

  ###Starting point###
  if (is.null(y0.start) && is.null(a.start)) {
    go <- 1
    while (go) {
      theta.start <- rnorm(2)*c(log.a.sd,y0.sd)+c(log.a.mean,y0.mean); theta.start[1] <- exp(theta.start[1])
      check <- pt.gen(x=theta.start[1]*(min.y-theta.start[2]), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      if (check >= min.prob) {
        go <- 0
      }
    }
    a <- theta.start[1]; phi <- log(a)
    y0 <- theta.start[2]
  } else {
    a <- a.start; phi <- log(a)
    y0 <- y0.start
  }
  prob.obs <- rep(1, n); prob.obs[ind.obs] <- pt.gen(x=a*(y[ind.obs]-y0), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  Psi <- U * (1 - as.numeric(ind.obs)/prob.obs)
  Psi.bar <- colMeans(Psi)
  Sigma <- Calc.EmpSigma(Psi = Psi, ind.use = (prob.obs >= min.prob))
  svd.Sigma <- svd(Sigma)
  
  for (i in 1:n.iter) {
    ##Metropolis step##
    #Update y0#
    cond.mean.y0 <- y0.mean + V[1,2]*(phi - log.a.mean)/V[1,1]
    y0.prop <- MH.y0(a = a, y0 = y0, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.y)
    ratio <- 0
    if (min.prob > 0) {
      tmp.quant <- min.y - qt(p=min.prob,df=t.df)/a
      ratio <-  -pnorm(q = tmp.quant, mean = y0.prop, sd = prop.y0.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = y0, sd = prop.y0.sd, log.p = T))
    }
    prob.obs.prop <- rep(1, n); prob.obs.prop[ind.obs] <- pt.gen(x=a*(y[ind.obs]-y0.prop), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    Psi.prop <- U * (1 - as.numeric(ind.obs)/prob.obs.prop)
    Psi.prop.bar <- colMeans(Psi.prop)
    Sigma.prop <- Calc.EmpSigma(Psi = Psi.prop, ind.use = (prob.obs.prop >= min.prob))
    svd.Sigma.prop <- svd(Sigma.prop)
    
    log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.bar)^2/svd.Sigma.prop$d ) - 1/2/cond.var.y0*(cond.mean.y0-y0.prop)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.bar)^2/svd.Sigma$d ) - 1/2/cond.var.y0*(cond.mean.y0-y0)^2 )
    if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
    if (log.post.diff + ratio > log(runif(1))) {  #Accept move
      y0 <- y0.prop
      prob.obs <- prob.obs.prop
      Psi <- Psi.prop
      Psi.bar <- Psi.prop.bar
      Sigma <- Sigma.prop
      svd.Sigma <- svd.Sigma.prop
    }
    
    #Update a#
    cond.mean.log.a <- log.a.mean + V[1,2]*(y0 - y0.mean)/V[2,2]
    a.prop <- MH.a(a = a, y0 = y0, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.y)
    phi.prop <- log(a.prop)
    ratio <- 0
    if (min.prob > 0 && min.y < y0 && min.prob < 1/2) {
      tmp.quant <- log(qt(p=min.prob,df=t.df)/(min.y-y0))
      ratio <-  -pnorm(q = tmp.quant, mean = phi.prop, sd = prop.a.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = phi, sd = prop.a.sd, log.p = T))
    }
    prob.obs.prop <- rep(1, n); prob.obs.prop[ind.obs] <- pt.gen(x=a.prop*(y[ind.obs]-y0), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    Psi.prop <- U * (1 - as.numeric(ind.obs)/prob.obs.prop)
    Psi.prop.bar <- colMeans(Psi.prop)
    Sigma.prop <- Calc.EmpSigma(Psi = Psi.prop, ind.use = (prob.obs.prop >= min.prob))
    svd.Sigma.prop <- svd(Sigma.prop)
    
    log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.bar)^2/svd.Sigma.prop$d ) - 1/2/cond.var.log.a*(cond.mean.log.a-phi.prop)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.bar)^2/svd.Sigma$d ) - 1/2/cond.var.log.a*(cond.mean.log.a-phi)^2 )
    if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
    if (log.post.diff + ratio > log(runif(1))) {  #Accept move
      a <- a.prop; phi <- log(a)
      prob.obs <- prob.obs.prop
      Psi <- Psi.prop
      Psi.bar <- Psi.prop.bar
      Sigma <- Sigma.prop
      svd.Sigma <- svd.Sigma.prop
    }
    
    ###Calculate statistics###
    if (save.every) {
      tmp.ind <- i/save.every
      if (floor(tmp.ind) == tmp.ind && tmp.ind <= n.save) {
        out$Samples$a[tmp.ind] <- a
        out$Samples$y0[tmp.ind] <- y0
      }
    }
    if (i > n.burn) {
      out$Post.Exp <- out$Post.Exp + 1/(n.iter-n.burn)*c(a, y0)
      Post.Exp2 <- Post.Exp2 + 1/(n.iter-n.burn)*cbind(c(a, y0))%*%rbind(c(a, y0))
      out$Post.Exp.W <- out$Post.Exp.W + 1/prob.obs/(n.iter-n.burn)
      out$Post.Exp.Pi <- out$Post.Exp.Pi + prob.obs/(n.iter-n.burn)
      Post.Exp.W2 <- Post.Exp.W2 + 1/prob.obs^2/(n.iter-n.burn)
    }
  }
  out$svd.Sigma <- svd.Sigma
  out$Psi.bar <- colMeans(Psi)
  out$Post.Var <- Post.Exp2 - cbind(out$Post.Exp)%*%rbind(out$Post.Exp) #Posteior variances of (a, y0)
  out$Post.Var.W <- Post.Exp.W2 - out$Post.Exp.W^2  #Posteior variances of pi^{-1}
  return(out)
}


###Calculate Sigma.hat###
Calc.EmpSigma <- function(Psi, ind.use=NULL) {
  if (!is.null(ind.use)) {
    Psi <- scale(Psi[ind.use,], center = T, scale = F)
  } else {
    Psi <- scale(Psi, center = T, scale = F)
  }
  return( 1/(nrow(Psi)-1) * t(Psi)%*%Psi )
}

#Estimate Sigma from the empirical distribution given missingness parameters#
Est.Sigma.Boot <- function(a, y0, y, U, ind.obs, t.df, n.Boot=1e3) {
  n <- length(y)
  return( Var.matrix( sapply(X = 1:n.Boot, function(b, a, y0, y, U, ind.obs, t.df, n){
    sample.ind.obs <- sample(x = ind.obs, size = n, replace = T, prob = 1/pt.gen(x = a*(y[ind.obs]-y0), df = t.df, p.min.1 = 0, p.min.2 = 0))
    y.try <- y[sample.ind.obs]
    U.try <- U[sample.ind.obs,]
    prob.try <- pt.gen(x = a*(y.try-y0), df = t.df, p.min.1 = 0, p.min.2 = 0)
    ind.obs.try <- rbinom(n = n, size = 1, prob = prob.try)
    return(sqrt(n) * colMeans(U.try * (1 - ind.obs.try/prob.try)))
  }, a=a, y0=y0, y=y, U=U, ind.obs=ind.obs, t.df=t.df, n=n) ) )
}

##Sample a with MH##
MH.a <- function(a, y0, t.df=4, p.min.1=0, p.min.2=0, met.sd, min.prob, min.y) {
  go <- 1
  while(go) {
    a.new <- a * exp(rnorm(1)*met.sd)
    if ( pt.gen(x=a.new*(min.y-y0), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2) >= min.prob ) {go <- 0}
  }
  return(a.new)
}

##Sample y0 with MH##
MH.y0 <- function(a, y0, t.df=4, p.min.1=0, p.min.2=0, met.sd, min.prob, min.y) {
  go <- 1
  while(go) {
    y0.new <- y0 + rnorm(1)*met.sd
    if ( pt.gen(x=a*(min.y-y0.new), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2) >= min.prob ) {go <- 0}
  }
  return(y0.new)
}



#######Old code########

Bayes.GMM_old <- function(Y, C, K.ind, p.min.1=0, p.min.2=0, t.df=4, y0.mean, log.a.mean, V, y0.start=NULL, a.start=NULL, prop.y0.sd=0.2, prop.a.sd=0.2, n.iter=1e4, n.burn=1e3, save.every=10, include.norm=T, min.prob=1/n) {
  max.store <- 1e4
  Y <- rbind(Y); C <- cbind(C); K.ind <- rbind(K.ind)
  n <- ncol(Y)
  p <- nrow(Y)
  K <- ncol(C)
  if (p.min.1 <= 0) {p.min.1 <- 0; p.min.2 <- 0}
  if (p.min.1 > 0 || min.prob <= 0) {min.prob <- 0}
  if (save.every <= 0) {
    save.every <- 0
  } else {
    save.every <- round(save.every)
    if (floor(n.iter/save.every) >= max.store) {save.every <- floor(n.iter/max.store)}
  }
  log.a.sd <- sqrt(V[1,1])
  y0.sd <- sqrt(V[2,2])
  cond.var.log.a <- V[1,1] - V[1,2]^2/V[2,2]
  cond.var.y0 <- V[2,2] - V[1,2]^2/V[1,1]
  
  #Define output#
  out <- list()
  out$Post.Exp <- matrix(0, nrow=p, ncol=2)  #Posterior expectation of (a, y0)
  Post.Exp2 <- lapply(X = 1:p, function(g){matrix(0,nrow=2,ncol=2)})  #Posterior expectation of (a, y0)^T(a, y0); a list of 2 x 2 matrices
  out$Post.Exp.W <- matrix(0, nrow=p, ncol=n)  #Posteior expectations of pi^{-1}
  out$Post.Exp.Pi <- matrix(0, nrow=p, ncol=n)  #Posteior expectations of pi
  Post.Exp.W2 <- matrix(0, nrow=p, ncol=n)  #Posteior expectations of pi^{-2}
  
  if (save.every) {
    out$Samples <- list(); n.save <- floor(n.iter/save.every); out$Samples$a <- matrix(NA, nrow=p, ncol=n.save); out$Samples$y0 <- matrix(NA, nrow=p, ncol=n.save)
  }
  A.vec <- rep(NA, p)
  Y0.vec <- rep(NA, p)
  
  for (i in 1:n.iter) {
    for (g in 1:p) {
      y.g <- Y[g,]
      ind.obs.g <- !is.na(y.g)
      min.yg <- min(y.g[ind.obs.g])
      
      ###Starting point###
      if (i == 1) {
        if (is.null(y0.start) && is.null(a.start)) {
          go <- 1
          while (go) {
            theta.start <- rnorm(2)*c(log.a.sd,y0.sd)+c(log.a.mean,y0.mean); theta.start[1] <- exp(theta.start[1])
            check <- pt.gen(x=theta.start[1]*(min.yg-theta.start[2]), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
            if (check >= min.prob) {
              go <- 0
            }
          }
          A.vec[g] <- theta.start[1]
          Y0.vec[g] <- theta.start[2]
        } else {
          A.vec[g] <- a.start[g]
          Y0.vec[g] <- y0.start[g]
        }
        if (save.every == 1) {out$Samples$a[g,i] <- A.vec[g]; out$Samples$y0[g,i] <- Y0.vec[g]}
        next
      }
      
      ###i > 1###
      a.g <- A.vec[g]; phi.g <- log(a.g); y0.g <- Y0.vec[g]
      U.g <- Get.U(C = C, K.ind = K.ind[g,])
      prob.obs.g <- rep(1, n); prob.obs.g[ind.obs.g] <- pt.gen(x=a.g*(y.g[ind.obs.g]-y0.g), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      Psi.g <- U.g * (1 - as.numeric(ind.obs.g)/prob.obs.g)
      Psi.g.bar <- colMeans(Psi.g)
      Sigma.g <- Calc.EmpSigma(Psi = Psi.g, ind.use = (prob.obs.g >= min.prob))
      svd.Sigma <- svd(Sigma.g)
      
      ##Metropolis step##
      #Update y0#
      cond.mean.y0 <- y0.mean + V[1,2]*(phi.g - log.a.mean)/V[1,1]
      y0.prop <- MH.y0(a = a.g, y0 = y0.g, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.yg)
      ratio.g <- 0
      if (min.prob > 0) {
        tmp.quant <- min.yg - qt(p=min.prob,df=t.df)/a.g
        ratio.g <-  -pnorm(q = tmp.quant, mean = y0.prop, sd = prop.y0.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = y0.g, sd = prop.y0.sd, log.p = T))
      }
      prob.obs.prop.g <- rep(1, n); prob.obs.prop.g[ind.obs.g] <- pt.gen(x=a.g*(y.g[ind.obs.g]-y0.prop), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      Psi.prop.g <- U.g * (1 - as.numeric(ind.obs.g)/prob.obs.prop.g)
      Psi.prop.g.bar <- colMeans(Psi.prop.g)
      Sigma.prop.g <- Calc.EmpSigma(Psi = Psi.prop.g, ind.use = (prob.obs.prop.g >= min.prob))
      svd.Sigma.prop <- svd(Sigma.prop.g)
      
      log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.g.bar)^2/svd.Sigma.prop$d ) - 1/2/cond.var.y0*(cond.mean.y0-y0.prop)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.g.bar)^2/svd.Sigma$d ) - 1/2/cond.var.y0*(cond.mean.y0-y0.g)^2 )
      if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
      if (log.post.diff + ratio.g > log(runif(1))) {  #Accept move
        y0.g <- y0.prop
        Y0.vec[g] <- y0.prop
        prob.obs.g <- prob.obs.prop.g
        Psi.g <- Psi.prop.g
        Psi.g.bar <- Psi.prop.g.bar
        Sigma.g <- Sigma.prop.g
        svd.Sigma <- svd.Sigma.prop
      }
      
      #Update a#
      cond.mean.log.a <- log.a.mean + V[1,2]*(y0.g - y0.mean)/V[1,1]
      a.prop <- MH.a(a = a.g, y0 = y0.g, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.yg)
      phi.prop <- log(a.prop)
      ratio.g <- 0
      if (min.prob > 0 && min.yg < y0.g) {
        tmp.quant <- log(qt(p=min.prob,df=t.df)/(min.yg-y0.g))
        ratio.g <-  -pnorm(q = tmp.quant, mean = phi.prop, sd = prop.a.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = phi.g, sd = prop.a.sd, log.p = T))
      }
      prob.obs.prop.g <- rep(1, n); prob.obs.prop.g[ind.obs.g] <- pt.gen(x=a.prop*(y.g[ind.obs.g]-y0.g), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      Psi.prop.g <- U.g * (1 - as.numeric(ind.obs.g)/prob.obs.prop.g)
      Psi.prop.g.bar <- colMeans(Psi.prop.g)
      Sigma.prop.g <- Calc.EmpSigma(Psi = Psi.prop.g, ind.use = (prob.obs.prop.g >= min.prob))
      svd.Sigma.prop <- svd(Sigma.prop.g)
      
      log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.g.bar)^2/svd.Sigma.prop$d ) - 1/2/cond.var.log.a*(cond.mean.log.a-phi.prop)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.g.bar)^2/svd.Sigma$d ) - 1/2/cond.var.log.a*(cond.mean.log.a-phi.g)^2 )
      if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
      if (log.post.diff + ratio.g > log(runif(1))) {  #Accept move
        A.vec[g] <- a.prop
        a.g <- a.prop
        prob.obs.g <- prob.obs.prop.g
      }
      
      ###Calculate statistics###
      if (save.every) {
        tmp.ind <- i/save.every
        if (floor(tmp.ind) == tmp.ind && tmp.ind <= n.save) {
          out$Samples$a[g,tmp.ind] <- a.g
          out$Samples$y0[g,tmp.ind] <- y0.g
        }
      }
      if (i > n.burn) {
        out$Post.Exp[g,] <- out$Post.Exp[g,] + 1/(n.iter-n.burn)*c(a.g, y0.g)
        Post.Exp2[[g]] <- Post.Exp2[[g]] + 1/(n.iter-n.burn)*cbind(c(a.g, y0.g))%*%rbind(c(a.g, y0.g))
        out$Post.Exp.W[g,] <- out$Post.Exp.W[g,] + 1/prob.obs.g/(n.iter-n.burn)
        out$Post.Exp.Pi[g,] <- out$Post.Exp.Pi[g,] + prob.obs.g/(n.iter-n.burn)
        Post.Exp.W2[g,] <- Post.Exp.W2[g,] + 1/prob.obs.g^2/(n.iter-n.burn)
      }
    }
  }
  
  out$Post.Var <- lapply(X = seq_len(p), function(g){Post.Exp2[[g]] - cbind(out$Post.Exp[g,])%*%rbind(out$Post.Exp[g,])})  #Posteior variances of (a, y0)
  out$Post.Var.W <- Post.Exp.W2 - out$Post.Exp.W^2  #Posteior variances of pi^{-1}
  return(out)
}


###Tests###
Test.Var <- function(X, K.ind, prob.obs, ind.obs, ind.use=NULL) {
  out <- list()
  U <- Get.U(C = cbind(X), K.ind = K.ind)
  Psi <- U; Psi[ind.obs,] <- Psi[ind.obs,] * (1 - 1/prob.obs)
  out$Psi.bar <- colMeans(Psi)
  out$Sigma <- Calc.EmpSigma(Psi = Psi, ind.use = ind.use)
  return(out)
}