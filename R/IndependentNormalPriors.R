#Bayesian GMM without knowing prior on missingness parameters
#This cannot be run in parallel
#This assumes the log(scale) and location parameters are independent N(mu,sigma^2), where mu is known and sigma^{-2} is Gamma(alpha, beta)

Bayes.GMM.ind <- function(Y, C, K.ind, p.min.1=0, p.min.2=0, t.df=4, y0.mean, log.a.mean, alpha.gamma = 1, beta.gamma = 1, y0.start=NULL, a.start=NULL, prop.y0.sd=0.2, prop.a.sd=0.2, n.iter=1e4, n.burn=1e3, save.every=10, ind.variance=NULL, include.norm=T, min.prob=1/n) {
  n <- ncol(Y)
  p <- nrow(Y)
  if (is.null(ind.variance)) {ind.variance <- 1:p}
  if (is.logical(ind.variance[1])) {ind.variance <- which(ind.variance == T)}
  p.var <- length(ind.variance)
  var.y0 <- 1
  var.log.a <- 1
  
  max.store <- 1e4
  if (p.min.1 <= 0) {p.min.1 <- 0; p.min.2 <- 0}
  if (p.min.1 > 0 || min.prob <= 0) {min.prob <- 0}
  if (save.every <= 0) {
    save.every <- 0
  } else {
    save.every <- round(save.every)
    if (floor(n.iter/save.every) >= max.store) {save.every <- floor(n.iter/max.store)}
  }
  
  
  if (is.null(y0.start)) {
    y0.start <- rnorm(p)*sqrt(var.y0) + y0.mean
  } else {
    var.y0 <- var.y0 <- 1/rgamma( n = 1, shape = p.var/2 + alpha.gamma, rate = 1/2*sum((y0.vec[ind.variance] - y0.mean)^2) + beta.gamma )
  }
  if (is.null(a.start)) {
    a.start <- exp(rnorm(p)*sqrt(var.log.a) + log.a.mean)
  } else {
    var.log.a <- 1/rgamma( n = 1, shape = p.var/2 + alpha.gamma, rate = 1/2*sum((log(a.vec[ind.variance]) - log.a.mean)^2) + beta.gamma )
  }
  
  #Initialize parameters#
  a.vec <- a.start
  y0.vec <- y0.start
  Prob.obs <- matrix(NA, nrow=p, ncol=n)
  
  #Initialize output#
  out <- list()
  if (save.every > 0) {
    out$a.saved <- matrix(NA, nrow=p, ncol=floor(n.iter/save.every))
    out$y0.saved <- matrix(NA, nrow=p, ncol=floor(n.iter/save.every))
    out$Var.saved <- matrix(NA, nrow=2, ncol=floor(n.iter/save.every))
  }
  out$Post.Exp <- matrix(0, nrow=p, ncol=2)
  Post.Exp2 <- lapply(X = 1:p, function(g){matrix(0,nrow=2,ncol=2)})
  out$Post.Exp.W <- matrix(0, nrow=p, ncol=n)
  out$Post.Exp.Pi <- matrix(0, nrow=p, ncol=n)
  Post.Exp.W2 <- matrix(0, nrow=p, ncol=n)
  out$Post.exp.Var <- c(0, 0)
  Post.exp.Var2 <- matrix(0, nrow=2, ncol=2)
  
  for (i in 1:n.iter) {
    for (g in 1:p) {
      a.g <- a.vec[g]; phi.g <- log(a.g)
      y0.g <- y0.vec[g]
      U.g <- Get.U(C = cbind(C), K.ind = K.ind[g,])
      y.g <- Y[g,]
      ind.obs.g <- !is.na(y.g)
      min.y.g <- min(y.g[ind.obs.g])
      
      #Get Psi, Psi.bar, Sigma and svd.Sigma
      if (is.na(Prob.obs[g,1])) {
        Prob.obs[g,] <- rep(1, n); Prob.obs[g,ind.obs.g] <- pt.gen(x=a.g*(y.g[ind.obs.g]-y0.g), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      }
      Psi <- U.g * (1 - as.numeric(ind.obs.g)/Prob.obs[g,])
      Psi.bar <- colMeans(Psi)
      Sigma <- Calc.EmpSigma(Psi = Psi, ind.use = (Prob.obs[g,] >= min.prob))
      svd.Sigma <- svd(Sigma)
      
      #y0#
      y0.prop.g <- MH.y0(a = a.g, y0 = y0.g, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.y.g)
      ratio <- 0
      if (min.prob > 0) {
        tmp.quant <- min.y.g - qt(p=min.prob,df=t.df)/a.g
        ratio <-  -pnorm(q = tmp.quant, mean = y0.prop.g, sd = prop.y0.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = y0.g, sd = prop.y0.sd, log.p = T))
      }
      prob.obs.prop <- rep(1, n); prob.obs.prop[ind.obs.g] <- pt.gen(x=a.g*(y.g[ind.obs.g]-y0.prop.g), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      Psi.prop <- U.g * (1 - as.numeric(ind.obs.g)/prob.obs.prop)
      Psi.prop.bar <- colMeans(Psi.prop)
      Sigma.prop <- Calc.EmpSigma(Psi = Psi.prop, ind.use = (prob.obs.prop >= min.prob))
      svd.Sigma.prop <- svd(Sigma.prop)
      
      log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.bar)^2/svd.Sigma.prop$d ) - 1/2/var.y0*(y0.mean-y0.prop.g)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.bar)^2/svd.Sigma$d ) - 1/2/var.y0*(y0.mean-y0.g)^2 )
      if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
      if (log.post.diff + ratio > log(runif(1))) {  #Accept move
        y0.vec[g] <- y0.prop.g
        y0.g <- y0.prop.g
        Prob.obs[g,] <- prob.obs.prop
        Psi <- Psi.prop
        Psi.bar <- Psi.prop.bar
        Sigma <- Sigma.prop
        svd.Sigma <- svd.Sigma.prop
      }
      
      #a#
      a.prop.g <- MH.a(a = a.g, y0 = y0.g, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, met.sd = prop.y0.sd, min.prob = min.prob, min.y = min.y.g)
      phi.prop.g <- log(a.prop.g)
      ratio <- 0
      if (min.prob > 0 && min.y.g < y0.g && min.prob < 1/2) {
        tmp.quant <- log(qt(p=min.prob,df=t.df)/(min.y.g-y0.g))
        ratio <-  -pnorm(q = tmp.quant, mean = phi.prop.g, sd = prop.a.sd, log.p = T) - (-pnorm(q = tmp.quant, mean = phi.g, sd = prop.a.sd, log.p = T))
      }
      prob.obs.prop <- rep(1, n); prob.obs.prop[ind.obs.g] <- pt.gen(x=a.prop.g*(y.g[ind.obs.g]-y0.g), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      Psi.prop <- U.g * (1 - as.numeric(ind.obs.g)/prob.obs.prop)
      Psi.prop.bar <- colMeans(Psi.prop)
      Sigma.prop <- Calc.EmpSigma(Psi = Psi.prop, ind.use = (prob.obs.prop >= min.prob))
      svd.Sigma.prop <- svd(Sigma.prop)
      
      log.post.diff <- ( -n/2*sum( as.vector(t(svd.Sigma.prop$u)%*%Psi.prop.bar)^2/svd.Sigma.prop$d ) - 1/2/var.log.a*(log.a.mean-phi.prop.g)^2 ) - ( -n/2*sum( as.vector(t(svd.Sigma$u)%*%Psi.bar)^2/svd.Sigma$d ) - 1/2/var.log.a*(log.a.mean-phi.g)^2 )
      if (include.norm) {log.post.diff <- log.post.diff - 1/2*sum(log(svd.Sigma.prop$d)) + 1/2*sum(log(svd.Sigma$d))}
      if (log.post.diff + ratio > log(runif(1))) {  #Accept move
        a.vec[g] <- a.prop.g
        Prob.obs[g,] <- prob.obs.prop
      }
      
      #Update estimators#
      if (i > n.burn) {
        out$Post.Exp[g,] <- out$Post.Exp[g,] + 1/(n.iter-n.burn)*c(a.vec[g], y0.vec[g])
        Post.Exp2[[g]] <- Post.Exp2[[g]] + 1/(n.iter-n.burn)*cbind(c(a.vec[g], y0.vec[g]))%*%rbind(c(a.vec[g], y0.vec[g]))
        out$Post.Exp.W[g,] <- out$Post.Exp.W[g,] + 1/Prob.obs[g,]/(n.iter-n.burn)
        out$Post.Exp.Pi[g,] <- out$Post.Exp.Pi[g,] + Prob.obs[g,]/(n.iter-n.burn)
        Post.Exp.W2[g,] <- Post.Exp.W2[g,] + 1/Prob.obs[g,]^2/(n.iter-n.burn)
      }
    }
    
    #Update var.y0 and var.log.a#
    var.y0 <- 1/rgamma( n = 1, shape = p.var/2 + alpha.gamma, rate = 1/2*sum((y0.vec[ind.variance] - y0.mean)^2) + beta.gamma )
    var.log.a <- 1/rgamma( n = 1, shape = p.var/2 + alpha.gamma, rate = 1/2*sum((log(a.vec[ind.variance]) - log.a.mean)^2) + beta.gamma )
    
    if (i > n.burn) {
      out$Post.exp.Var <- out$Post.exp.Var + 1/(n.iter-n.burn)*c(var.log.a, var.y0)
      Post.exp.Var2 <- Post.exp.Var2 + 1/(n.iter-n.burn)*cbind(c(var.log.a, var.y0))%*%rbind(c(var.log.a, var.y0))
    }
    
    if (save.every > 0) {
      if (floor(i/save.every) == i/save.every && floor(i/save.every) <= ncol(out$a.saved)) {
        out$a.saved[,floor(i/save.every)] <- a.vec
        out$y0.saved[,floor(i/save.every)] <- y0.vec
        out$Var.saved[,floor(i/save.every)] <- c(var.log.a, var.y0)
      }
    }
  }
  out$Post.Var <- lapply(X = 1:p, function(g){Post.Exp2[[g]] - cbind(out$Post.Exp[g,])%*%rbind(out$Post.Exp[g,])})
  out$Post.Var.W <- Post.Exp.W2 - out$Post.Exp.W^2
  out$Post.Var.Var <- Post.exp.Var2 - cbind(out$Post.exp.Var)%*%rbind(out$Post.exp.Var)
}