library(sva)
library(qvalue)
#####Do confounder correctin with missing data#####
#Input is data matrix, X, Z, K and missingness mechanism estimated with EstimateMissing

#' Estimate latent covariates and coefficients of interest
#'
#' Estimate latent covariates and estimate/do inference on the coefficients of interest in a multivariate linear model with stabilized inverse probability weighting (sIPW) using estimated missingness mechansim.
#'
#' @param Y a \code{p} x \code{n} data matrix of log2-transformed metabolite intensities, where \code{p} = #of metabolites and \code{n} = #of samples. Missing values should be left as \code{NA}.
#' @param X a \code{n} x \code{d} matrix of covariates of interest (i.e. disease status).
#' @param Z a \code{n} x \code{r} matrix of observed nuisance covariates (i.e. the intercept, observed technical factors, etc.)
#' @param K The number of latent covariates (i.e. C is a \code{n} x \code{K} matrix). If unspecified, it is estimated using sva::num.sv applied to the metabolites with complete data.
#' @param Miss.Mech The missingness mechansim object returned by \code{EstimateMissing}.
#' @param ind.samples A logical or numeric vector of samples to be considered in the analysis. For example, if disease status were only measured in a subset of the patients, this would be the samples with a recorded disease status. Default is no missing samples.
#' @param est.Beta A logical indicating whether or not to estimate/do inference on coefficients of interest. If \code{F}, only the latent covariates are estimated. Default, and recommended value, is \code{T}.
#' @param max.miss.perp Maximum fraction of missing data a metabolite is allowed to have to be used to calculate the part of C perpendicular to X. Defaults to 0.5.
#' @param max.miss.image Maximum fraction of missing data a metabolite is allowed to have to be used to calculate the part of C in the image of X. Defaults to 0.5.
#' 
#' @return A list \item{C.iter}{The estimate of the \code{n} x \code{K} matrix of latent covariates.} \item{Beta.iter}{The estimate of the \code{p} x \code{d} matrix of coefficients of interest.} \item{p.t.iter}{The \code{p} x \code{d} matrix of p-values for the coefficients of interest.} \item{Var.beta.iter}{A length \code{p} list of the \code{d} x \code{d} estimates for Var(Beta.iter)} \item{t.iter}{A \code{p} x \code{d} matrix of t-statistics, defined as Beta.iter/SE(Beta.iter)} \item{p.f.iter}{A length \code{p} vector of F-statistic p-values for the null hypothesis that X has no effect on metabolite intensity. It is only returned if \code{d} > 1.} \item{L}{The estimate of the \code{p} x \code{K} matrix of coefficients for the latent covariates.} \item{Beta.naive}{The estimate of the \code{p} x \code{d} matrix of coefficients of interest that ignores C. This should ONLY be used for comparison.} \item{p.t.naive}{The \code{p} x \code{d} matrix of p-values for the coefficients of interest that ignore C. This should ONLY be used for comparison.}
#'
#' @export
CC.Missing <- function(Y, X, Z=NULL, K=NULL, Miss.Mech, ind.samples=NULL, max.miss.perp=0.5, max.miss.image=0.5, BH.min=NULL, method = c("sIPW", "IPW"), include.updates=T, est.Beta=T, return.nuis=F, refine.C=F, p.refine.both=F, return.all=F, return.mu=F) {
  method <- match.arg(method, c("sIPW", "IPW"))
  max.miss.C <- Miss.Mech$max.miss.C
  X <- cbind(X); Z <- cbind(Z)
  out <- list()
  
  #Fix observed samples#
  if (is.null(ind.samples)) {
    ind.samples <- 1:nrow(X)
  } else {
    if (is.logical(ind.samples[1])) {
      length.samples <- sum(ind.samples)
    } else {
      length.samples <- length(ind.samples)
    }
    if (ncol(Y) > length.samples) {Y <- Y[,ind.samples]}
    if (nrow(X) > length.samples) {X <- cbind(X[ind.samples,])}
    if (!is.null(Z)) {
      if (nrow(Z) > length.samples) {Z <- cbind(Z[ind.samples,])}
    }
  }
  
  d <- ncol(X)
  p <- nrow(Y)
  n <- ncol(Y)
  Prob.Missing <- rowMeans(is.na(Y))
  out$X <- X; out$Z <- Z; out$method <- method
  ind.miss.all <- !is.na(Miss.Mech$Theta.Miss[,1]) & !is.na(Miss.Mech$Theta.Miss[,2])  #Metabolites with missing data
  if (is.null(BH.min)) {BH.min <- Miss.Mech$BH.min}
  tmp <- BH.proc(p = Miss.Mech$Pvalue.value, alpha = BH.min)

  out$flag <- !(tmp == FALSE & !is.na(tmp)) & Prob.Missing > max.miss.C  #All flagged metabolites
  out$flag.analyzed <- !(tmp == FALSE & !is.na(tmp)) & ind.miss.all    #Metabolites that have reported p-values, but the missingness mechanism is suspect
  
  ####Initial estimate for C.perp####
  #If K is null, estimate it with sva::num.sv
  if (is.null(K)) {K <- sva::num.sv(dat = Y[Prob.Missing==0,], mod = cbind(out$X,out$Z))}; out$K <- K
  if (include.updates) {cat(paste("Estimating", K, "latent factors...", collapse = ""))}
  out.C <- EstC.0(Y = Y, K = K, Cov = cbind(X,Z), max.miss = max.miss.C, max.iter = 800, n.repeat.Sigma = 1)
  out$C.perp <- out.C$C
  if (include.updates) {cat("done\n")}
  
  #Estimate C without including missing metabolites#
  tmp.cov.ignore <- cbind(out$X,out$C.perp,out$Z)
  Y1L.tmp <- t(apply(X = Y[Prob.Missing <= max.miss.C,], MARGIN = 1, function(y){ind.obs.y<-!is.na(y); (solve(t(tmp.cov.ignore[ind.obs.y,])%*%tmp.cov.ignore[ind.obs.y,],t(tmp.cov.ignore[ind.obs.y,])%*%y[ind.obs.y]))[1:(d+K)]}))
  Omega.ignore.tmp <- t(cbind(Y1L.tmp[,1:d]))%*%Y1L.tmp[,(d+1):(d+K)]%*%solve(t(Y1L.tmp[,(d+1):(d+K)])%*%Y1L.tmp[,(d+1):(d+K)])
  if (d == 1) {Omega.ignore.tmp <- rbind(as.vector(Omega.ignore.tmp))}
  if (return.all) {out$C.ignore.miss <- X%*%Omega.ignore.tmp + out$C.perp}
  
  ####Determine Weights####
  if (method == "sIPW") {
    Weights <- Miss.Mech$Post.W[,ind.samples] * Miss.Mech$Pi.MAR[,ind.samples]
    Var.Weights <- Miss.Mech$Pi.MAR[,ind.samples]^2*( Miss.Mech$Post.W[,ind.samples]^2 + Miss.Mech$Post.VarW[,ind.samples] )
  }
  if (method == "IPW") {
    Weights <- Miss.Mech$Post.W[,ind.samples]
    Var.Weights <- Miss.Mech$Post.W[,ind.samples]^2 + Miss.Mech$Post.VarW[,ind.samples]
  }
  Weights[is.na(Weights)] <- 1
  
  ####Refine estimate for C.perp####
  out$C.perp.old <- NULL
  if (max.miss.perp > max.miss.C) {
    if (include.updates) {cat(paste("Refining estimate for latent factors with missingness mech...", collapse = ""))}
    ind.perp <- Prob.Missing <= max.miss.C | (Miss.Mech$Ind.Confident & Prob.Missing <= max.miss.perp)
    out.perp <- Est.Cperp.Weights(Y = Y, Cov = cbind(X,Z), C.start = out$C.perp, Weights = Weights, ind.use.miss = ind.perp, max.miss.C = max.miss.C)
    out$C.perp.old <- out$C.perp
    out$n.iter.refine <- out.perp$n.iter
    out$C.perp <- out.perp$C.perp
    if (include.updates) {cat("done\n")}
  } else {
    ind.perp <- Prob.Missing <= max.miss.C
  }
  if (!est.Beta) {return(out)}
  
  ####Estimate Omega####
  if (include.updates) {cat(paste("Performing estimation and inference...", collapse = ""))}
  Cov <- cbind(X,out$C.perp,Z)
  BLG.naive <- matrix(NA, nrow=p, ncol=ncol(Cov))
  Var.naive <- matrix(NA, nrow=p, ncol=d)
  Var.C <- vector(mode = "list", length = p)
  
  #Complete data
  ind.complete <- Prob.Missing <= max.miss.C
  tmp <- lapply(X = which(ind.complete==T), function(g){ ind.obs <- !is.na(Y[g,])
            y <- Y[g,ind.obs]
            cov.tmp <- Cov[ind.obs,]
            hess <- solve(t(cov.tmp)%*%cov.tmp)
            tt <- list(); tt$coef <- hess%*%t(cov.tmp)%*%y
            sigma2 <- 1/(length(y)-ncol(cov.tmp))*sum((y-cov.tmp%*%tt$coef)^2)
            tt$var.C <- sigma2*hess[(d+1):(d+K),(d+1):(d+K)]
            tt$var.x <- sigma2*diag(hess)[1:d]
            return(tt) })
  BLG.naive[ind.complete,] <- t(sapply(tmp,function(x){x$coef}))
  Var.naive[ind.complete,] <- t(sapply(tmp,function(x){x$var.x}))
  Var.C[ind.complete] <- lapply(tmp,function(x){x$var.C})
  
  #Missing data
  ind.nmar <- Prob.Missing > max.miss.C & ind.miss.all
  if (sum(ind.nmar) > 0) {
    tmp <- lapply(X = which(ind.nmar==T), function(g){ ind.obs <- !is.na(Y[g,])
                y <- Y[g,ind.obs]
                w <- Weights[g,ind.obs]
                var.w <- Var.Weights[g,ind.obs]
                cov.tmp <- Cov[ind.obs,]
                hess <- solve(t(cov.tmp*w)%*%cov.tmp)
                tt <- list(); tt$coef <- hess%*%t(cov.tmp*w)%*%y
                Score <- cov.tmp*(y-as.vector(cov.tmp%*%tt$coef))
                Score.vc <- t(sapply(X = 1:nrow(cov.tmp), function(i){ as.vector(solve( diag(1,ncol(cov.tmp),ncol(cov.tmp)) - w[i]*cbind(cov.tmp[i,])%*%rbind(cov.tmp[i,])%*%hess, Score[i,] )) }))
                tmp.var <- hess%*%( t(Score.vc*var.w)%*%Score.vc )%*%hess
                tt$var.C <- tmp.var[(d+1):(d+K),(d+1):(d+K)]
                tt$var.x <- diag(tmp.var)[1:d]
                return(tt) })
    BLG.naive[ind.nmar,] <- t(sapply(tmp,function(x){x$coef}))
    Var.naive[ind.nmar,] <- t(sapply(tmp,function(x){x$var.x}))
    Var.C[ind.nmar] <- lapply(tmp,function(x){x$var.C})
  }
  out$L <- BLG.naive[,(d+1):(d+K)]
  
  #Estimate Omega#
  if (max.miss.image > max.miss.C) {
    ind.image <- Prob.Missing <= max.miss.C | (Miss.Mech$Ind.Confident & Prob.Missing <= max.miss.image) 
  } else {
    ind.image <- Prob.Missing <= max.miss.C
  }
  B.naive.Omega <- cbind(BLG.naive[ind.image,1:d])
  L.Omega <- cbind(BLG.naive[ind.image,(d+1):(d+K)])
  Var.Omega <- cbind(Var.naive[ind.image,])
  Var.C.Omega <- Var.C[ind.image]
  
  ####Estimate L'Sigma.invL####
  out$LtSL <- 1/nrow(L.Omega) * t(L.Omega/Var.Omega[,1])%*%L.Omega
  out$RtSR <- 1/nrow(L.Omega) * Reduce("+", lapply(X = 1:sum(ind.image), function(g){1/Var.Omega[g,1]*Var.C.Omega[[g]]}))

  out$Omega <- matrix(NA, nrow=K, ncol=d)
  out$Omega.iter <- matrix(NA, nrow=K, ncol=d)
  out$Omega.corr <- matrix(NA, nrow=K, ncol=d)
  out$Omega.pvalue <- rep(NA, d)
  out$Omega.pvalue.corr <- rep(NA, d)
  if (!is.null(Z)) {
    X.tilde <- cbind(X - Z%*%solve(t(Z)%*%Z,t(Z)%*%X))
  }
  for (r in 1:d) {
    LSL.r <- t(L.Omega/Var.Omega[,r])%*%L.Omega
    Omega.r <- solve(LSL.r, t(L.Omega/Var.Omega[,r])%*%B.naive.Omega[,r])
    Omega.r.corr <- solve(LSL.r - Reduce("+", lapply(X = 1:sum(ind.image), function(g){1/Var.Omega[g,r]*Var.C.Omega[[g]]})), t(L.Omega/Var.Omega[,r])%*%B.naive.Omega[,r])
    out$Omega[,r] <- Omega.r
    out$Omega.corr[,r] <- Omega.r.corr
    out$Omega.iter[,r] <- Iterate.Omega(Y = Y, Y1 = B.naive.Omega[,r] / sqrt(Var.Omega[,r]), L = cbind(L.Omega / sqrt(Var.Omega[,r])), X = out$X[,r], Z = cbind(out$X[,-r],out$Z), Omega = Omega.r, C.perp = out$C.perp, B.iter = 3, q.thresh = 0.1, Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C, ind.image = ind.image)
    if (is.null(Z)) {
      out$Omega.pvalue[r] <- pchisq(sum(Omega.r^2)*sum(X[,r]^2), df = K, lower.tail = F)
      out$Omega.pvalue.corr[r] <- pchisq(sum(Omega.r.corr^2)*sum(X[,r]^2), df = K, lower.tail = F)
    } else {
      out$Omega.pvalue[r] <- pchisq(sum(Omega.r^2)*sum(X.tilde[,r]^2), df = K, lower.tail = F)
      out$Omega.pvalue.corr[r] <- pchisq(sum(Omega.r.corr^2)*sum(X.tilde[,r]^2), df = K, lower.tail = F)
    }
  }
  out$C <- out$X %*% t(out$Omega) + out$C.perp
  out$C.corr <- out$X %*% t(out$Omega.corr) + out$C.perp
  out$C.iter <- out$X %*% t(out$Omega.iter) + out$C.perp
  
  if (refine.C && return.all) {  #Refine C using IRW-SVA-like procedure with out$C.corr as a starting point
    tmp.IRW <- IRW.C(Y = Y, X = out$X, C = out$C.corr, Z = out$Z, Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.perp = ind.perp, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C, p.naive = F)
    out$C.IRW <- tmp.IRW$C
    out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$C.IRW,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C)
    out$p.t.IRW <- out.tmp$p.t
    out$Beta.IRW <- out.tmp$Beta
    out$Var.beta.IRW <- out.tmp$Var.beta
    out$Sigma.IRW <- out.tmp$Sigma
    
    if (d == 1 && p.refine.both) {
      tmp.IRW <- IRW.C(Y = Y, X = out$X, C = out$C.corr, Z = out$Z, Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.perp = ind.perp, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C, p.naive = T)
      out$C.IRW.naive <- tmp.IRW$C
      out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$C.IRW.naive,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C)
      out$p.t.IRW.naive <- out.tmp$p.t
      out$Beta.IRW.naive <- out.tmp$Beta
      out$Var.beta.IRW.naive <- out.tmp$Var.beta
      out$Sigma.IRW.naive <- out.tmp$Sigma
    }
  }
  
  ####Estimate and inference on B without iterating####
  if (return.all) {
    out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$C,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C)
    out$t <- out.tmp$t
    out$p.t <- out.tmp$p.t
    out$f <- out.tmp$f
    out$p.f <- out.tmp$p.f
    out$Beta <- out.tmp$Beta
    out$Var.beta <- out.tmp$Var.beta
  }
  
  ####Iterative estimate for Omega####
  out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$C.iter,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C, return.mu=return.mu, return.nuis=return.nuis)
  out$t.iter <- out.tmp$t
  out$p.t.iter <- out.tmp$p.t
  out$f.iter <- out.tmp$f
  out$p.f.iter <- out.tmp$p.f
  if (return.nuis) {
    out$Beta.iter <- out.tmp$Beta[,1:d]
    out$L.iter <- out.tmp$Beta[,(d+1):(d+K)]
    if (!is.null(out$Z)) {out$Nuis.iter <- out.tmp$Beta[,(d+K+1):ncol(cbind(out$X,out$C.iter,out$Z))]}
  } else {
    out$Beta.iter <- out.tmp$Beta
  }
  out$Var.beta.iter <- out.tmp$Var.beta  
  out$Sigma.iter <- out.tmp$Sigma
  if (return.mu) {out$Mu <- out.tmp$Mu}
  
  ######Inflated estimate for Omega######
  if (return.all) {
    out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$C.corr,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C)
    out$p.t.corr <- out.tmp$p.t
    out$Beta.corr <- out.tmp$Beta
    out$Var.beta.corr <- out.tmp$Var.beta
    out$Sigma.corr <- out.tmp$Sigma
  }
  
  ######Ignore C######
  out.tmp.naive <- Estimation.IPW(Y = Y, Cov = cbind(out$X,out$Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:d, ind.analyze = ind.miss.all | Prob.Missing <= max.miss.C)
  out$p.t.naive <- out.tmp.naive$p.t.naive  #Ignore C
  out$Beta.naive <- out.tmp.naive$Beta.naive
  
  ######Ignore Missingness######
  if (ncol(X) == 1 && return.all) {
    out$p.t.ignore.miss <- rep(NA, nrow(Y))
    out$Beta.ignore.miss <- rep(NA, nrow(Y))
    out$Var.beta.ignore.miss <- rep(NA, nrow(Y))
    out$Sigma.ignore.miss <- rep(NA, nrow(Y))
    tmp.ignore <- Ignore.Missing(Y = Y[ind.miss.all | Prob.Missing <= max.miss.C,], Cov = cbind(out$X,out$C.ignore.miss,out$Z))
    out$p.t.ignore.miss[ind.miss.all | Prob.Missing <= max.miss.C] <- tmp.ignore$p
    out$Beta.ignore.miss[ind.miss.all | Prob.Missing <= max.miss.C] <- tmp.ignore$beta
    out$Var.beta.ignore.miss[ind.miss.all | Prob.Missing <= max.miss.C] <- tmp.ignore$sd.beta^2
    out$Sigma.ignore.miss[ind.miss.all | Prob.Missing <= max.miss.C] <- tmp.ignore$Sigma
  }

  
  if (include.updates) {cat("done\n")}
  return(out)
}



###Estimate C.perp with missing data weights###

Est.Cperp.Weights <- function(Y, Cov=NULL, C.start, Weights, ind.use.miss, max.miss.C=0.05, max.iter=800, tol=1e-6) {
  Frac.missing <- apply(Y, 1, function(y) {mean(is.na(y))})
  ind.use <- Frac.missing <= max.miss.C | ind.use.miss
  Frac.missing <- Frac.missing[ind.use]
  Y <- Y[ind.use,]
  Weights <- Weights[ind.use,]
  K <- ncol(C.start)
  if (!is.null(Cov)) {
    Cov <- cbind(Cov)
    d <- ncol(Cov)
    Cov.perp <- diag(nrow(Cov)) - Cov%*%solve(t(Cov)%*%Cov,t(Cov))
  } else {
    d <- 0
  }
  
  #Start function after pruning#
  n <- ncol(Y)
  p <- nrow(Y)
  
  C.perp <- C.start
  Cov.total <- cbind(Cov, C.perp)
  Ind.Obs <- !is.na(Y)
  log.like <- -1e16
  log.like.vec <- rep(NA, max.iter+1)
  for (i in 1:max.iter) {
    BL <- t(sapply(X = 1:p, function(g){ ind.obs <- !is.na(Y[g,]); y <- Y[g,ind.obs]; w <- Weights[g,ind.obs]; solve(t(Cov.total[ind.obs,]*w)%*%Cov.total[ind.obs,],t(Cov.total[ind.obs,]*w)%*%y) }))
    L <- cbind(BL[,(d+1):(d+K)])
    
    if (!is.null(Cov)) {
      B <- cbind(BL[,1:d])
      C.perp <- t(sapply(X = 1:n, function(i){ ind.obs <- !is.na(Y[,i]); y <- Y[ind.obs,i]-as.vector(B[ind.obs,]%*%cbind(Cov[i,])); w <- Weights[ind.obs,i]; solve(t(L[ind.obs,]*w)%*%L[ind.obs,],t(L[ind.obs,]*w)%*%y) }))
    } else {
      C.perp <- t(sapply(X = 1:n, function(i){ ind.obs <- !is.na(Y[,i]); y <- Y[ind.obs,i]; w <- Weights[ind.obs,i]; solve(t(L[ind.obs,]*w)%*%L[ind.obs,],t(L[ind.obs,]*w)%*%y) }))
    }
    Cov.total <- cbind(Cov, C.perp)
    
    ##Compute log-likelihood##
    Mu <- BL %*% t(Cov.total)
    log.like.new <- -sum(Weights[Ind.Obs] * (Y[Ind.Obs]-Mu[Ind.Obs])^2)
    log.like.vec[i] <- log.like.new
    if (abs(log.like.new - log.like)/abs(log.like) <= tol) {
      break
    }
    log.like <- log.like.new
  }
  
  if (!is.null(Cov)) {
    C.perp <- sqrt(n-d)*svd(Cov.perp %*% C.perp)$u
  } else {
    C.perp <- sqrt(n)*svd(C.perp)$u
  }
  
  return(list(C.perp=C.perp, n.iter=i, log.like=log.like.vec[1:i]))
}


######Estimation with IPW######

Estimation.IPW <- function(Y=Y, Cov=Cov, Weights, Var.Weights, max.miss.C=0.05, ind.int=NULL, ind.analyze, return.mu=F, return.nuis=F) {
  n <- ncol(Y)
  p <- nrow(Y)
  r <- ncol(Cov)
  if (is.null(ind.int) || return.nuis) {ind.int <- 1:r}
  out <- list()
  
  Prob.Missing <- apply(X = Y, MARGIN = 1, function(x) {mean(is.na(x))})
  Beta <- matrix(NA, nrow=p, ncol=r)
  Beta.naive <- matrix(NA, nrow=p, ncol=r)
  Var.list <- vector(mode = "list", length = p)   #Array of variance matrices
  Var.naive.list <- vector(mode = "list", length = p)
  out$Sigma <- rep(NA, p)
  if (return.mu) {out$Mu <- matrix(NA, nrow=p, ncol=n)}
  
  #Missing at random#
  ind.mar <- Prob.Missing <= max.miss.C
  tmp.mar <- lapply(X = which(ind.mar==T), function(g){ ind.obs <- !is.na(Y[g,])
          y <- Y[g,ind.obs]
          cov.tmp <- cbind(Cov[ind.obs,])
          hess <- solve(t(cov.tmp)%*%cov.tmp)
          tt <- list(); tt$coef <- hess%*%t(cov.tmp)%*%y
          tt$sigma2 <- 1/(length(y)-ncol(cov.tmp))*sum((y-cov.tmp%*%tt$coef)^2)
          tt$var <- tt$sigma2*hess
  return(tt) })
  Beta[ind.mar,] <- t(sapply(tmp.mar, function(x){x$coef}))
  Beta.naive[ind.mar,] <- Beta[ind.mar,]
  Var.list[ind.mar] <- lapply(tmp.mar, function(x){x$var})
  Var.naive.list[ind.mar] <- Var.list[ind.mar]
  out$Sigma[ind.mar] <- unlist(lapply(tmp.mar, function(x){x$sigma2}))
  if (return.mu) {out$Mu[ind.mar,] <- Beta[ind.mar,]%*%t(Cov)}
  
  #Missing not at random#
  ind.nmar <- Prob.Missing > max.miss.C & ind.analyze
  Q.cov <- Compute.Orthog.Proj(Cov)
  if (sum(ind.nmar) > 0) {
    tmp.nmar <- lapply(X = which(ind.nmar==T), function(g){ ind.obs <- !is.na(Y[g,])
                  y <- Y[g,ind.obs]
                  w <- Weights[g,ind.obs]
                  var.w <- Var.Weights[g,ind.obs]
                  cov.tmp <- cbind(Cov[ind.obs,])
                  hess <- solve(t(cov.tmp*w)%*%cov.tmp); hess.naive <- solve(t(cov.tmp)%*%cov.tmp)
                  tt <- list(); tt$coef <- hess%*%t(cov.tmp*w)%*%y; tt$coef.naive <- hess.naive%*%t(cov.tmp)%*%y
                  tt$mu <- 
                  resids.ipw <- y-as.vector(cov.tmp%*%tt$coef)
                  Score <- cov.tmp*resids.ipw
                  Score.naive <- cov.tmp*(y-as.vector(cov.tmp%*%tt$coef.naive))
                  
                  Score.vc <- t(sapply(X = 1:nrow(cov.tmp), function(i){ as.vector(solve( diag(1,ncol(cov.tmp),ncol(cov.tmp)) - w[i]*cbind(cov.tmp[i,])%*%rbind(cov.tmp[i,])%*%hess, Score[i,] )) }))
                  
                  tt$var <- hess%*%( t(Score*var.w)%*%Score )%*%hess; tt$var.vc <- hess%*%( t(Score.vc*var.w)%*%Score.vc )%*%hess; tt$var.naive <- sum((y-cov.tmp%*%tt$coef.naive)^2)/(nrow(cov.tmp)-ncol(cov.tmp))*hess.naive
                  
                  #tt$sigma2 <- sum(w*(y-as.vector(cov.tmp%*%tt$coef))^2)/( sum(w) - sum(diag( hess%*%(t(cov.tmp*w^2)%*%cov.tmp) )) )  #This gives better estimates of sigma2 in practice
                  hat.ipw <- apply(X = cov.tmp*sqrt(w), MARGIN = 1, function(x,A){sum(x*as.vector(A%*%x))}, A=hess)
                  tt$sigma2 <- sum( w*resids.ipw^2 / (1-hat.ipw)^2 )/sum(w)
                  return(tt) })
    Beta[ind.nmar,] <- t(sapply(tmp.nmar,function(x){x$coef}))
    Beta.naive[ind.nmar,] <- t(sapply(tmp.nmar,function(x){x$coef.naive}))
    Var.list[ind.nmar] <- lapply(tmp.nmar,function(x){x$var.vc})
    Var.naive.list[ind.nmar] <- lapply(tmp.nmar,function(x){x$var.naive})
    out$Sigma[ind.nmar] <- unlist(lapply(tmp.nmar, function(x){x$sigma2}))
    if (return.mu) {out$Mu[ind.nmar,] <- Beta[ind.nmar,]%*%t(Cov)}
  }
  
  
  #t#
  out$t <- matrix(NA, nrow=p, ncol=length(ind.int))
  out$t[ind.analyze,] <- t(rbind(sapply( X = which(ind.analyze==T), function(g){ Beta[g,ind.int] / sqrt(diag(Var.list[[g]])[ind.int]) } )))
  out$p.t <- 2*pnorm(-abs(out$t))
  
  out$p.t.naive <- matrix(NA, nrow=p, ncol=length(ind.int))
  out$p.t.naive[ind.analyze,] <- 2*pnorm(-abs( t(rbind(sapply( X = which(ind.analyze==T), function(g){ Beta.naive[g,ind.int] / sqrt(diag(Var.naive.list[[g]])[ind.int]) } ))) ))

  if (length(ind.int) > 1) {
    out$f <- rep(NA, nrow=p)
    out$f[ind.analyze] <- unlist(lapply(X = which(ind.analyze==T), function(g){ sum(Beta[g,ind.int] * solve(Var.list[[g]][ind.int,ind.int], Beta[g,ind.int])) }))/length(ind.int)
    out$p.f <- pchisq(q = length(ind.int)*out$f, df = length(ind.int), lower.tail = FALSE)
  } else {
    out$f <- NULL
    out$p.f <- NULL
  }
  out$Beta <- Beta[,ind.int]
  out$Beta.naive <- Beta.naive[,ind.int]
  out$Var.beta <- lapply(Var.list, function(x){if(is.null(x)){return(NA)}; if(is.na(c(x)[1])){return(NA)}; x[ind.int,ind.int]})
  return(out)
}


Ignore.Missing <- function(Y, Cov, ind.int=1) {
  z <- lapply(X = 1:nrow(Y), function(g){ out <- list(); ind.id <- !is.na(Y[g,]); y <- Y[g,ind.id]; cov <- cbind(Cov[ind.id,]); hess.g <- solve(t(cov)%*%cov); out$beta <- as.vector(hess.g%*%t(cov)%*%y); sigma2 <- sum((y-cov%*%out$beta)^2)/(nrow(cov)-ncol(cov)); out$sigma2 <- sigma2; out$z <- out$beta/sqrt(sigma2*diag(hess.g)); return(out) })
  return(list(p=2*pnorm(-abs(unlist(lapply(z, function(x){x$z[ind.int]})))), beta=unlist(lapply(z, function(x){x$beta[ind.int]})), Sigma=unlist(lapply(z, function(x){x$sigma2})), sd.beta=unlist(lapply(z, function(x){x$beta[ind.int]/x$z[ind.int]}))))
}

Compute.Orthog.Proj <- function(X) {
  X <- cbind(X)
  qr.X <- qr(X)
  Q <- qr.Q(qr.X, complete = T)[,(qr.X$rank+1):nrow(X)]
  return(Q%*%t(Q))
}

###Iterative estimation of Omega###
Iterate.Omega <- function(Y, Y1, L, X, Z=NULL, Omega, C.perp, B.iter=3, q.thresh=0.1, Weights, Var.Weights, max.miss.C, ind.analyze, ind.image) {
  X <- cbind(X)
  Omega <- cbind(Omega)
  for (j in 1:B.iter) {
    C <- C.perp + X%*%t(Omega)
    tmp.j <- Estimation.IPW(Y = Y, Cov = cbind(X,C,Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:ncol(X), ind.analyze = ind.analyze)
    q.j <- qvalue::qvalue(as.vector(tmp.j$p.t)[ind.image])
    ind.j <- q.j$qvalues > q.thresh & !is.na(q.j$qvalues)
    if (sum(ind.j) == sum(ind.image)) {return(as.vector(Omega))}
    Omega <- cbind(solve(t(L[ind.j,]) %*% L[ind.j,], t(L[ind.j,]) %*% Y1[ind.j]))
  }
  return(as.vector(Omega))
}


###Iteratively-reweighted BCconf###
#X contains the covariates of interest. For now, I will assume dim(X) = 1
#Z contains the nuisance covariates, like the intercept
#C.start contains a starting point for the latent factors
#The output is an estimate for C, which is orthogonal to Z
IRW.C <- function(Y, X, Z, C, Weights, Var.Weights, max.miss.C, ind.perp, ind.analyze, B.iter=5, p.naive=F) {  
  X <- cbind(X)
  Z <- cbind(Z)
  C <- cbind(C)
  
  d <- ncol(X)
  K <- ncol(C)
  p <- nrow(Y)
  
  w.C <- rep(NA, p)
  for (i in 1:B.iter) {
    out.tmp <- Estimation.IPW(Y = Y, Cov = cbind(X,C,Z), Weights = Weights, Var.Weights = Var.Weights, max.miss.C = max.miss.C, ind.int = 1:(d+K), ind.analyze = ind.analyze)
    
    #Pvalues for B#
    if (d == 1) {
      if (p.naive) {
        p.X <- out.tmp$p.t.naive[,1]
      } else {
        p.X <- out.tmp$p.t[,1]
      }
    } else {
      p.X <- pchisq(q = unlist(lapply( 1:p, function(g){ if (is.na(out.tmp$Beta[g,1])){return(NA)}; sum(out.tmp$Beta[g,1:d]*as.vector(solve(out.tmp$Var.beta[[g]][1:d,1:d],out.tmp$Beta[g,1:d]))) } )), df = d, lower.tail = F)
    }
    
    #Pvalues for L#
    if (K > 1) {
      p.C <- pchisq(q = unlist(lapply( 1:p, function(g){ if (is.na(out.tmp$Beta[g,1])){return(NA)}; sum(out.tmp$Beta[g,(d+1):(K+d)]*as.vector(solve(out.tmp$Var.beta[[g]][(d+1):(K+d),(d+1):(K+d)],out.tmp$Beta[g,(d+1):(K+d)]))) } )), df = K, lower.tail = F)
    } else {
      p.C <- 2*pnorm(-abs( out.tmp$Beta[,d+K] / sqrt( unlist(lapply(out.tmp$Var.beta,function(x){if (is.na(c(x)[1])){return(NA)}; x[d+K,d+K]})) ) ))
    }
    
    #P(B = 0, L != 0 | Data)#
    w.C[ind.analyze] <- edge.lfdr(p.X[ind.analyze])*(1-edge.lfdr(p.C[ind.analyze]))  #Factor analysis weights
    
    #Initial estimate for C using only fully-observed metabolites#
    out.C <- EstC.0(Y = Y * w.C, K = K, Cov = Z, max.miss = max.miss.C, max.iter = 800, n.repeat.Sigma = 1)
    
    #Refine estimate for C using only fully-observed metabolites#
    out.C <- Est.Cperp.Weights(Y = Y * w.C, Cov = Z, C.start = out.C$C, Weights = Weights, ind.use.miss = ind.perp, max.miss.C = max.miss.C)
    C <- out.C$C.perp
    if (i < B.iter) {w.C <- rep(NA, p)}
  }
  return(list(C=C, w.C=w.C))
}
  

####Compute P(beta = 0 | Data) from sva-devel/R/helper.R on github####

edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  n <- length(p)
  transf <- match.arg(transf)
  
  if(transf=="probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1-eps)
    x <- qnorm(p)
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x <- log((p+eps)/(1-p+eps))
    myd <- density(x, adjust=adj)
    mys <- smooth.spline(x=myd$x, y=myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1+exp(x))^2
    lfdr <- pi0 * dx/y
  }
  
  if(trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if(monotone) {	
    lfdr <- lfdr[order(p)]
    lfdr <- mono(lfdr)
    lfdr <- lfdr[rank(p)]
  }
  
  return(lfdr)
}

mono <- function(lfdr){
  .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
}