#####Estimate T_4 CDF missingness mechanism#####

Opt.GMM.t.grid <- function(A.grid, Pos.Grid, y, C=NULL, K.ind=NULL, U=NULL, W=NULL, t.df=4, n.iter=1, return.all=F, p.min.1=0, p.min.2=0) {
  t.df <- t.df[1]
  if (is.null(U)) {U <- cbind(Get.U(C, K.ind=K.ind))}
  ind.obs <- !is.na(y)
  y.orig <- y
  y <- y[ind.obs]
  n <- length(y)
  
  m.grid <- length(Pos.Grid)
  if (is.null(W)) {W <- 1/n*t(U)%*%U}
  Out.grid <- t(sapply(A.grid, function(a){
    prob.obs <- pt.gen(cbind(rep(1,m.grid))%*%rbind(a*y) - cbind(a*Pos.Grid)%*%rbind(rep(1,n)), df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2) #m.grad x n matrix

    G <- cbind(rep(1,m.grid))%*%rbind(colSums(U)) - (1/prob.obs)%*%U[ind.obs,]  #m.grid x k
    out.a <- 1/2*rowSums((G%*%W)*G)
    out.a[is.nan(out.a)] <- Inf
    return(out.a)
  }))  #A.grid x Pos.Grid
  ind.a <- which.min(apply(X = Out.grid, MARGIN = 1, min))
  ind.y <- which.min(apply(X = Out.grid, MARGIN = 2, min))
  a <- A.grid[ind.a]
  y0 <- Pos.Grid[ind.y]*A.grid[ind.a]
  if (n.iter == 1) {
    if (ind.a == 1 || ind.a == length(A.grid) || ind.y == 1 || ind.y == length(Pos.Grid)) {return(list(a=a, y0=y0, out.value=min(Out.grid), Out.grid=Out.grid))}
    #out <- optim(par = c(a,y0), fn = func.GMM.norm, gr = Grad.GMM.norm, method = "BFGS", y=y.orig, U=U, W=W)
    out <- constrOptim(theta = c(a,y0), f = func.GMM.t, grad = Grad.GMM.t, ui = rbind(c(1,0), c(-1,0), c(0,1), c(0,-1)), ci = c(A.grid[ind.a-1], -A.grid[ind.a+1], Pos.Grid[ind.y-1]*A.grid[ind.a-1], -Pos.Grid[ind.y+1]*A.grid[ind.a+1]), y=y.orig, U=U, W=W, t.df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    if (return.all) {
      return(list(a=out$par[1], y0=out$par[2], out.value=out$value, W=W, a.init=a, y0.init=y0, out.init=Out.grid[ind.a,ind.y], Out.grid=Out.grid, A.grid=A.grid, Pos.Grid=Pos.Grid))
    }
    return(list(a=out$par[1], y0=out$par[2], out.value=out$value, W=W, a.init=a, y0.init=y0, out.init=Out.grid[ind.a,ind.y]))
  }
  
  G <- U; G[ind.obs,] <- U[ind.obs,] * (1 - 1/pt.gen(a*y-y0, df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2))
  return(Opt.GMM.t.grid(A.grid=A.grid, Pos.Grid=Pos.Grid, y=y.orig, W=solve(1/length(y.orig)*t(G)%*%G), U = U, n.iter=n.iter-1, t.df=t.df, return.all=return.all, p.min.1 = p.min.1, p.min.2 = p.min.2))
}

##Compute value at optimum from Opt.GMM.t.grid##

Val.Opt.t <- function(a, y0, y, C, K.ind = 1:2, t.df=4, p.min.1=0, p.min.2=0) {
  U <- cbind(Get.U(C, K.ind=K.ind))
  ind.obs <- !is.na(y)
  G <- U; G[ind.obs,] <- U[ind.obs,] * (1 - 1/pt.gen(x = a*(y[ind.obs]-y0), df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2))
  G.bar <- colMeans(G)
  n <- nrow(U)
  return( n*sum(G.bar*as.vector(solve(Var.matrix(X = G),G.bar))) )
}

#####Functions#####

#theta is a vector

func.GMM.t <- function(theta, y, U, W, t.df=4, p.min.1 = 0, p.min.2 = 0) {
  U <- cbind(U)
  a <- theta[1]
  y0 <- theta[2]
  ind.obs <- !is.na(y)
  y <- y[ind.obs]
  
  prob.obs <- pt.gen(a*y-y0, df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  prob.obs[prob.obs <= 0] <- 1e-8
  if (ncol(U) == 1) {
    g <- sum(U) - sum(U[ind.obs,]/prob.obs)
    return(1/2 * g^2)
  } else {
    g <- colSums(U) - colSums(U[ind.obs,]/prob.obs)
    return(1/2 * sum( g*(W %*% g) ))
  }
}


#####Gradients#####

#Gradient for first and second iteration#
Grad.GMM.t<- function(theta, y, U, W, t.df=4, p.min.1 = 0, p.min.2 = 0) {
  U <- cbind(U)
  a <- theta[1]
  y0 <- theta[2]
  n <- length(y)
  ind.obs <- !is.na(y)  
  y <- y[ind.obs]
  prob.obs <- pt.gen(a*y-y0, df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  prob.obs[prob.obs <= 0] <- 1e-8
  
  Grad.pi <- cbind(y*dt.gen(a*y-y0, df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2), -dt.gen(a*y-y0, df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2))
  Grad.g <- t(Grad.pi/prob.obs^2) %*% U[ind.obs,]  #2 x k
  if (ncol(U) == 1) {
    return( as.vector(Grad.g * (sum(U) - sum(U[ind.obs,]/prob.obs))) )
  } else {
    return( as.vector(Grad.g %*% W %*% (colSums(U) - colSums(U[ind.obs,]/prob.obs))) )
  }
}



#######Variance of estimator#######
#Computes variance when we parametrize the function as Phi{a(y - mu)}, where a*mu=y0
Var.GMM.t <- function(theta, y, C, K.ind, W=NULL, t.df=4, p.min.1 = 0, p.min.2 = 0) {  #According to numericDeriv, this is correct
  t.df <- t.df[1]
  if (is.null(W)) {W <- diag(length(K.ind)+1)}
  U <- Get.U(C = C, K.ind = K.ind)
  ind.obs <- !is.na(y)
  y.obs <- y[ind.obs]
  prob.obs <- pt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  Psi <- U; Psi[ind.obs,] <- Psi[ind.obs,]*(1-1/prob.obs)     #k-vector
  V <- Calc.EmpSigma(Psi = Psi) * (nrow(Psi)-1)
  Gamma <- t(U[ind.obs,]/prob.obs^2*dt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)) %*% cbind(y.obs-theta[2]/theta[1], rep(-theta[1],sum(ind.obs)))  #A k x p matrix
  
  Middle <- t(Gamma) %*% W %*% V %*% W %*% Gamma
  H1 <- t(Gamma) %*% W %*% Gamma
  
  return(list( Asy.Var=solve(H1,Middle)%*%solve(H1), H1=H1, Asy.Var.2step=solve(t(Gamma) %*% solve(V) %*% Gamma) ))  #Total Hessian = H1+H2 when the gradient of f is 0, where H2 should be negligible as n -> infinity, which appears to be the case.
}


Var.InvPi2 <- function(theta, y, C, K.ind, W=NULL, t.df=4, transform=T, p.min.1 = 0, p.min.2 = 0) {
  if (is.null(W)) {W <- diag(length(K.ind)+1)}
  if (transform) {theta[2] <- theta[2]*theta[1]}
  V <- Var.GMM.t(theta = theta, y = y, C = C, K.ind = K.ind, W = W, t.df = t.df)$Asy.Var
  ind.obs <- !is.na(y)
  y.obs <- y[ind.obs]
  n.obs <- sum(ind.obs)
  y.stand <- theta[1]*y.obs - theta[2]
  #prob.obs <- pnorm(y.stand); min.y <- min(y.stand[prob.obs >= 0.2]); y.stand[y.stand < min.y] <- min.y
  prob.obs <- pt.gen(x = y.stand, df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  Grad.pi <- cbind(y.obs, rep(-1,n.obs)) * dt.gen(x = y.stand, df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  S <- -y.stand*dt.gen(x = y.stand, df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  Out <- rep(NA, n.obs)
  for (j in 1:n.obs) {
    tmp.j <- (3/prob.obs[j]^4*cbind(Grad.pi[j,])%*%rbind(Grad.pi[j,]) - S[j]/prob.obs[j]^3*matrix(c(y.obs[j]^2,-y.obs[j],-y.obs[j],1), nrow=2, ncol=2)) %*% V
    Out[j] <- sum(diag(tmp.j))
  }
  return(Out)
}


pt.gen <- function(x, df=4, p.min.1=0, p.min.2=0) {
  if (p.min.1 == 0) {
    return(pt(x, df=df))
  }
  p.min.2 <- max(p.min.2, p.min.1)
  q.1 <- qt(p.min.1, df=df)
  q.2 <- qt(p.min.2, df=df)
  tmp <- pt(x, df=df)
  tmp[tmp <= p.min.1] <- p.min.1
  if (p.min.1 == p.min.2) {return(tmp)}
  #Linear interpolation#
  #tmp[tmp >= p.min.1 & tmp <= p.min.2] <- (p.min.2 - p.min.1)/(q.2 - q.1) * (x[tmp >= p.min.1 & tmp <= p.min.2] - q.1) + p.min.1
  #tmp[tmp < p.min.1] <- p.min.1
  #return(tmp)
  #Spline interpolation#
  mm <- solve(cbind(c((q.2-q.1)^3,3*(q.2-q.1)^2),c((q.2-q.1)^2,2*(q.2-q.1))),c(p.min.2-p.min.1,dt(q.2,df=df)))
  tmp[tmp > p.min.1 & tmp < p.min.2] <- mm[1]*(x[tmp > p.min.1 & tmp < p.min.2] - q.1)^3 + mm[2]*(x[tmp > p.min.1 & tmp < p.min.2] - q.1)^2 + p.min.1
  return(tmp)
}

dt.gen <- function(x, df=4, p.min.1=0, p.min.2=0) {
  if (p.min.1 == 0) {
    return(dt(x, df=df))
  }
  p.min.2 <- max(p.min.2, p.min.1)
  q.1 <- qt(p.min.1, df=df)
  q.2 <- qt(p.min.2, df=df)
  tmp <- dt(x, df=df)
  tmp[x <= q.1] <- 0
  if (p.min.1 == p.min.2) {return(tmp)}
  mm <- solve(cbind(c((q.2-q.1)^3,3*(q.2-q.1)^2),c((q.2-q.1)^2,2*(q.2-q.1))),c(p.min.2-p.min.1,dt(q.2,df=df)))
  tmp[x > q.1 & x < q.2] <- 3*mm[1]*(x[x > q.1 & x < q.2] - q.1)^2 + 2*mm[2]*(x[x > q.1 & x < q.2] - q.1)
  #tmp[x >= q.1 & x <= q.2] <- (p.min.2 - p.min.1)/(q.2 - q.1)
  return(tmp)
}


####Empirical Bayes to estimate hyperparameters####

Emp.Bayes.MuSigma <- function(Mu.g, Var.g, middle=c("mean", "median"), shift.var=0, refine.mu=T, simple.average=F) {
  middle <- match.arg(middle, c("mean", "median"))
  ind.use <- !is.na(Mu.g) & !is.na(Var.g)
  Mu.g <- Mu.g[ind.use]
  Var.g <- Var.g[ind.use]
  if (middle == "median") {
    mu <- median(Mu.g)
  } else {
    if (simple.average){
      mu <- mean(Mu.g)
    } else {
      mu <- sum(Mu.g/(Var.g+shift.var))/sum(1/(Var.g+shift.var))
    }
  }
  
  Mu.g <- Mu.g - mu
  Phi <- 10^(seq(-2,4,by=0.05))
  like <- unlist(lapply(X = Phi, FUN = function(phi){ sum(Mu.g^2/Var.g/(1+Var.g*phi) - log(1+Var.g*phi) + log(phi)) }))
  like.inf <- -sum(log(Var.g))
  ind.max <- which.max(like)
  
  #Refine estimate for phi#
  if (like.inf < max(like) && ind.max > 1) {
    Phi2 <- seq(Phi[ind.max-1], Phi[ind.max+1], by=min(0.01,(Phi[ind.max+1]-Phi[ind.max-1])/5))
    like2 <- unlist(lapply(X = Phi2, FUN = function(phi){ sum(Mu.g^2/Var.g/(1+Var.g*phi) - log(1+Var.g*phi) + log(phi)) }))
    if (refine.mu && middle == "mean" && !simple.average) {
      return(Emp.Bayes.MuSigma(Mu.g = Mu.g+mu, Var.g = Var.g, middle = "mean", shift.var = 1/Phi2[which.max(like2)], refine.mu = F))
    }
    return(list(mu=mu, sd=1/sqrt(Phi2[which.max(like2)]), like=like, Phi=Phi))
  }
  if (like.inf >= max(like)) {
    return(list(mu=mu, sd=0, like=like, Phi=Phi))
  }
  if (ind.max == 1) {
    return(list(mu=mu, sd=1/sqrt(Phi[1]), like=like, Phi=Phi))
  }
}

##EB on both a and y0 simultaneously##

Emp.Bayes.MuSigma.Both <- function(Mu, Var, simple.average=F) {
  out <- list()
  ind.use <- unlist(lapply(X = 1:nrow(Mu), function(g){if (sum(is.na(Mu[g,])) > 0){return(F)}; if (is.na(Var[[g]])){return(F)}; if (class( try(expr = {solve(Var[[g]])}, silent = T) ) == "try-error"){return(F)}; tmp.s <- svd(Var[[g]]); if (tmp.s$d[1]/tmp.s$d[2] > 1e6){return(F)}; return(T)}))
  Mu <- Mu[ind.use,]
  out$ind.use <- ind.use
  Var <- lapply(which(ind.use==T), function(g){Var[[g]]})
  p <- nrow(Mu)
  
  Emp.a <- Emp.Bayes.MuSigma(Mu.g = Mu[,1], Var.g = unlist(lapply(1:p,function(g){Var[[g]][1,1]})), middle = "mean", refine.mu = T, simple.average = simple.average)
  Emp.y0 <- Emp.Bayes.MuSigma(Mu.g = Mu[,2], Var.g = unlist(lapply(1:p,function(g){Var[[g]][2,2]})), middle = "mean", refine.mu = T, simple.average = simple.average)
  mu <- c(Emp.a$mu, Emp.y0$mu)
  out$mu <- mu
  var.a <- max(Emp.a$sd^2, 1e-2)
  var.y0 <- max(Emp.y0$sd^2, 1e-2)
  
  out.EB <- optim(par = c(1/var.a,0,1/var.y0), fn = f.EmpBayes, gr = grad.EmpBayes, hessian = T, Mu = Mu-cbind(rep(1,p))%*%rbind(mu), Var = Var)
  D <- matrix(NA, nrow=2, ncol=2)
  D[1,1] <- out.EB$par[1]; D[2,2] <- out.EB$par[3]; D[1,2] <- out.EB$par[2]; D[2,1] <- out.EB$par[2]
  out$V <- solve(D)
  out$Hessian <- out.EB$hessian
  return(out)
}

f.EmpBayes <- function(theta, Mu, Var) {
  D <- matrix(NA, nrow=2, ncol=2)
  D[1,1] <- theta[1]; D[2,2] <- theta[3]; D[1,2] <- theta[2]; D[2,1] <- theta[2]
  p <- nrow(Mu)
  like <- p*sum(log(eigen(D,symmetric = T)$values)) - sum(unlist(lapply(1:p,function(g){sum(log(svd(solve(Var[[g]])+D)$d))}))) + sum(unlist(lapply(1:p,function(g){ tmp.V <- Var[[g]]+Var[[g]]%*%D%*%Var[[g]]; tmp.svd <- svd(tmp.V); mu.rotate <- as.vector(t(tmp.svd$u)%*%Mu[g,]); sum(mu.rotate^2/tmp.svd$d) })))
  return(-like)
}

grad.EmpBayes <- function(theta, Mu, Var) {
  D <- matrix(NA, nrow=2, ncol=2)
  D[1,1] <- theta[1]; D[2,2] <- theta[3]; D[1,2] <- theta[2]; D[2,1] <- theta[2]
  p <- nrow(Mu)
  grad <- p*solve(D) - Reduce( "+", lapply(1:p,function(g){solve(D+solve(Var[[g]]))}) ) + Reduce( "+", lapply(1:p,function(g){ tmp.v <- solve(D+solve(Var[[g]]),solve(Var[[g]],Mu[g,])); tmp.v%*%t(tmp.v) }) )
  return(-c(grad[1,1], grad[1,2], grad[2,2]))
}



#####Boostrapped estimator for the variance#####

Boot.Var <- function(y, C=NULL, K.ind=1:2, U=NULL, a, y0, n.iter=2, t.df=4, p.min.1=0, p.min.2=0, n.Boot=100, y0.plus=2, y0.minus=2) {
  t.df <- t.df[1]
  y0.max <- 35
  y0.min <- 5
  n <- length(y)
  if (is.null(U)) {U <- cbind(Get.U(C, K.ind=K.ind))}
  A.grid <- seq(max(a/5, 0.1),min(5*a,7),by=0.05)
  Pos.Grid <- seq(max(y0-y0.minus,y0.min),min(y0+y0.plus,y0.max),by=0.25)
  
  cl <- makeCluster(max(detectCores() - 1, 1))
  clusterExport(cl, c("U", "A.grid", "Pos.Grid", "y", "n.iter", "t.df", "p.min.1", "p.min.2", "n"), envir=environment())
  clusterEvalQ(cl = cl, expr = source("../R/EstimateMissingnessGMM_t.R"))
  out.parallel <- parSapply(cl = cl, X = 1:n.Boot, function(b) {
    sample.b <- sample(x = 1:n, size = n, replace = T)
    out.b <- Opt.GMM.t.grid(A.grid = A.grid, Pos.Grid = Pos.Grid, y = y[sample.b], U = U[sample.b,], n.iter = n.iter, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    return(c(out.b$a, out.b$y0/out.b$a))
  })
  stopCluster(cl)

  return(list(a.boot=out.parallel[1,], y0.boot=out.parallel[2,], mean.a=mean(out.parallel[1,]), mean.y0=mean(out.parallel[2,]), mean.loga=mean(log(out.parallel[1,])), Var=Var.matrix(out.parallel), Var.log=Var.matrix(rbind(log(out.parallel[1,]),out.parallel[2,]))))
}

##Empirical Likelihood Boostrapped estimator for J-statistic##
Boot.J <- function(y, C, theta, K.ind=1:2, value.opt=NULL, t.df=4, p.min.1=0, p.min.2=0, n.boot=500, tol.EL=1e-8, y0.plus=2, y0.minus=2, plotit=F, alpha.wolfe=1/2, c.wolfe=0.1) {
  a <- theta[1]
  y0 <- theta[2]
  if (is.null(value.opt)) {value.opt <- Val.Opt.t(a = a, y0 = y0, y = y, C = C, K.ind = K.ind, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)}
  y0.max <- 35
  y0.min <- 5
  n <- length(y)
  ind.obs <- !is.na(y)
  U <- Get.U(C = C, K.ind = K.ind)
  Psi <- U; Psi[ind.obs,] <- Psi[ind.obs,] * (1 - 1/pt.gen(x = a*(y[ind.obs]-y0), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2))
  
  #Estimate EL weights#
  lambda <- rep(0,ncol(Psi))
  max.iter.EL <- 1e4
  for (i in 1:max.iter.EL) {
    weights <- 1/(1 - as.vector(Psi %*% lambda))  #EL weights
    grad.i <- as.vector( t(Psi)%*%weights )
    if (max(abs(grad.i)) < tol.EL) {
      break
    }
    hess.i <- t(Psi * weights^2) %*% Psi
    dir.i <-  as.vector(solve(hess.i,grad.i))
    f.i <- sum(log(weights))
    for (j in 1:20) {
      lambda.j <- lambda - dir.i
      weights.j <- 1/(1 - as.vector(Psi %*% lambda.j))
      if (min(weights.j) <= 0) {
        dir.i <- alpha.wolfe*dir.i
        next
      }
      if (sum(log(weights.j)) <= f.i - c.wolfe*sum(dir.i*grad.i)) {
        lambda <- lambda.j
        break
      } else {
        dir.i <- alpha.wolfe*dir.i
      }
    }
  }
  
  #Boostrap with estimated weights#
  A.grid <- seq(max(a/5, 0.1),min(5*a,7),by=0.05)
  Pos.Grid <- seq(max(y0-y0.minus,y0.min),min(y0+y0.plus,y0.max),by=0.25)
  value <- rep(NA, n.boot)
  for (b in 1:n.boot) {
    sample.b <- sample(x = 1:n, size = n, replace = T, prob = weights)
    y.b <- y[sample.b]
    C.b <- C[sample.b,] 
    out.boot <- Opt.GMM.t.grid(A.grid = A.grid, Pos.Grid = Pos.Grid, y = y.b, C = C.b, K.ind = K.ind, n.iter = 2, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    value[b] <- Val.Opt.t(a = out.boot$a, y0 = out.boot$y0/out.boot$a, y = y.b, C = C.b, K.ind = K.ind, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    #value[b] <- 2*out.boot$out.value/n
  }
  if (sum(value>=value.opt) < 5) {
    beta.start <- mean(value)/var(value)
    alpha.start <- beta.start*mean(value)
    out.gamma <- constrOptim(theta = c(alpha.start,beta.start), f = mlog.like.Gamma, grad = mgrad.like.Gamma, ui = diag(2), ci = c(1e-2,1e-2), mu.log=mean(log(value)), mu.x=mean(value))
    if (plotit) {
      hist(value, prob=T)
      lines(seq(0,max(value),length=1000),dgamma(seq(0,max(value),length=1000),out.gamma$par[1],out.gamma$par[2]),col="red")
    }
    return(list(value.opt=value.opt,value.boot=value,p.value.boot=(sum(value>=value.opt)+1)/(n.boot+1),p.value=pgamma(q = value.opt, shape = out.gamma$par[1], rate = out.gamma$par[2], lower.tail = F)))
  } else {
    return(list(value.opt=value.opt, value.boot=value, p.value=(sum(value>=value.opt)+1)/(n.boot+1)))
  }
  
}

mlog.like.Gamma <- function(theta, mu.log, mu.x) {
  alpha <- theta[1]; beta <- theta[2]
  out <- alpha*log(beta) - lgamma(alpha) + (alpha-1)*mu.log - beta*mu.x
  return(-out)
}

mgrad.like.Gamma <- function(theta, mu.log, mu.x) {
  alpha <- theta[1]; beta <- theta[2]
  grad.alpha <- log(beta) - digamma(alpha) + mu.log
  grad.beta <- alpha/beta - mu.x
  return(-c(grad.alpha,grad.beta))
}

Var.matrix <- function(X) {
  if (nrow(X) > ncol(X)) {
    X <- t(X)
  }
  n <- ncol(X)
  mu <- cbind(rowMeans(X))
  return( 1/(n-1)*X%*%t(X) - n/(n-1)*mu%*%t(mu) )
}

##U function##

Get.U <- function(C, K.ind=c(1,2)) {
  cbind(C)
  if (is.null(K.ind)) {return(cbind(rep(1,nrow(C))))}
  cbind(rep(1,nrow(C)), C[,K.ind])
}




############################ OLD FUNCTIONS ############################

Var.GMM.t.old <- function(theta, y, C, K.ind, W=NULL, t.df=4, p.min.1 = 0, p.min.2 = 0) {  #According to numericDeriv, this is correct
  t.df <- t.df[1]
  if (is.null(W)) {W <- diag(length(K.ind)+1)}
  U <- Get.U(C = C, K.ind = K.ind)
  ind.obs <- !is.na(y)
  y.obs <- y[ind.obs]
  prob.obs <- pt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  Psi <- U; Psi[ind.obs,] <- Psi[ind.obs,]*(1-1/prob.obs)     #k-vector
  V <- t(Psi)%*%Psi   #k x k
  Dh <- matrix(0, nrow=2, ncol=2); Dh[1,1] <- 1; Dh[1,2] <- theta[2]/theta[1]; Dh[2,2] <- theta[1]
  Gamma <- t(U[ind.obs,]/prob.obs^2*dt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)) %*% cbind(y.obs, rep(-1,sum(ind.obs))) %*% t(Dh)   #A k x p matrix
  Middle <- t(Gamma) %*% W %*% V %*% W %*% Gamma
  
  H1 <- t(Gamma) %*% W %*% Gamma
  Const <- as.vector(U[ind.obs,] %*% W %*% colSums(Psi))
  S <- -(theta[1]*y.obs-theta[2])*dt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  tmp.mat <- cbind(y.obs, rep(-1,sum(ind.obs)))
  H2 <- -2*t(tmp.mat / prob.obs^3 * Const * dt.gen(theta[1]*y.obs-theta[2], df=t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)^2) %*% tmp.mat
  H2[1,1] <- H2[1,1] + sum(S*y.obs^2/prob.obs^2)
  H2[1,2] <- H2[1,2] - sum(S*y.obs/prob.obs^2)
  H2[2,1] <- H2[2,1] - sum(S*y.obs/prob.obs^2)
  H2[2,2] <- H2[2,2] + sum(S/prob.obs^2)
  
  return(list( Asy.Var=solve(H1,Middle)%*%solve(H1), H1=H1, H2=Dh%*%H2%*%t(Dh) ))  #Total Hessian = H1+H2 when the gradient of f is 0, where H2 should be negligible as n -> infinity, which appears to be the case.
}
