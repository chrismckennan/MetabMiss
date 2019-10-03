library(parallel)
library(sva)
####Estimate missingness mechanism from metabolite data matrix####

#' Estimate the metabolite-dependent missingness mechanisms
#'
#' Estimate the metabolite-dependent missingness mechanisms with a hierarchical generalized method of moments (GMM). This function only has to be run once per metabolite dataset and the output should be stored with the metabolite data. The user need only specify \code{Y} and maybe \code{K}, although the default \code{K = 10} should suffice.
#' 
#' @param Y a \code{p} x \code{n} data matrix of log2-transformed metabolite intensities, where \code{p} = #of metabolites and \code{n} = #of samples. Missing values should be left as \code{NA}.
#' @param K a number >= 2. This gives the number of latent covariates to use to estimate the missingness mechanism. We recommend using \code{Num.Instruments} to estimate it. The default is 10. \code{K} and \code{Y} are the only variables that must be specified.
#' @param n_cores The number of cores to use. The default is the number of maximum number of usable cores - 1.
#' @param max.missing.consider The maximum fraction of missing data a metabolite is allowed to have. Missingness mechanisms will NOT be estimated for metabolites with more missing data than this. Default, and recommended value, is 0.5
#' @param max.miss.C Maximum fraction of missing data a metabolite can have to ignore the missingness mechanism in downstream estimation and inference. The default, and recommended value, is 0.05.
#' @param max.iter.C Maximum number of iterations to estimate the latent covariates C. Default is 400 and should not be changed.
#' @param t.df The missingness mechanism is the CDF of a scaled and cetered T-distribution with t.df degrees of freedom. The default, and recommended value, is 4
#' @param n.boot.J The number of bootstrap samples to compute the J-statistics. The defualt is 150.
#' @param Cov An optional n x d matrix of covariates. It is recommended the user not specify anything other than the intercept. The default is the intercept.
#' @param n.K.GMM Number of additional terms (besides the intercept) to be considered in GMM when estimating the missingness mechanism. The default, and recommendend value, is 2. If changed, this must be >= 2
#' @param Model.Pvalue A logical value. If \code{T}, a missingness model P-value is computed. The default, and recommended value, is \code{T}.
#'
#' @return A list that should be save immediately. It can be used directly as input into CC.Missing to estimate latent factors and the coefficients of interest in a multivariate linear model. \item{Post.Theta}{\code{p} x 2 matrix containing the posterior expectations of the missingness scale and location parameters (a,y0) for each metabolite. Returns \code{NA} for metabolites without a missingness mechansim.} \item{Post.Var}{A list of \code{p} 2x2 matrices containing the posterior variances for (a,y0) for each metabolite.} \item{Post.W}{A \code{p}x\code{n} containing the posterior expectations of 1/P(Metab is observed | y, a, y0), where the expectation is taken with respect to (a,y0) | y} \item{Post.VarW}{A \code{p}x\code{n} containing the posterior variances of 1/P(Metab is observed | y, a, y0), where the expectation is taken with respect to (a,y0) | y} \item{Post.Pi}{A \code{p}x\code{n} containing the posterior expectations of P(Metab is observed | y, a, y0), where the expectation is taken with respect to (a,y0) | y} \item{Pi.MAR}{A \code{p}x\code{n} containing estimate of P(Metab is observed | Latent covariates). This helps stabilize the inverse probability weights in downstream estimation.} \item{Theta.Miss}{\code{p} x 2 matrix with the estimates of the unshrunk GMM scale and location parameters a, y0 for each metabolite's missingness mechanism. If a missingness mechansism was not estimated, returns \code{NA}.} \item{Pvalue.value}{The J-test P-value that tests the null hypothesis H_0: Missingness mechanism is correct} \item{Ind.Confident}{A logical \code{p}-vector containing the indices of metabolites whose missingness mechanisms we are confident in.} \item{Emp.Bayes.loga}{Empirical Bayes estimate of E(log(a))} \item{Emp.Bayes.y0}{Empirical Bayes estimate of E(y0)}
#'
#' @export
EstimateMissing <- function(Y, K=10, max.missing.consider=0.5, Cov = NULL, max.miss.C = 0.05, n_cores=NULL, max.iter.C = 400, n.repeat.Sigma.C = 1, n.K.GMM = 2, min.a=0.1, max.a=7, min.y0=10, max.y0=30, t.df=4, p.min.1=0, p.min.2=0, n.boot.J=150, Model.Pvalue=T, BH.analyze.min=0.2, min.quant.5=5, shrink.Est=T, prop.y0.sd = 0.2, prop.a.sd = 0.2, n.iter.MCMC = 2e4, n.burn.MCMC = 1e3, min.prob.MCMC = 1/n, Bayes.est = c("EmpBayes", "FullBayes", "FullBayes_ind"), simple.average.EB=F) {
  p <- nrow(Y)
  n <- ncol(Y)
  Prob.Missing <- apply(X = Y, MARGIN = 1, function(x) {mean(is.na(x))})
  if (is.null(Cov)) {
    Cov <- cbind(rep(1,n))
  } else {Cov <- cbind(Cov)}
  out <- list()
  out$max.miss.C <- max.miss.C
  
  ###Estimate C###
  cat(paste("Estimating", K, "latent factors...", collapse = ""))
  out.C <- EstC.0(Y = Y, K = K, Cov = Cov, max.miss = max.miss.C, max.iter = max.iter.C, n.repeat.Sigma = n.repeat.Sigma.C)
  C <- out.C$C
  W <- out.C$W
  out$C <- out.C$C
  cat("done\n")  
  
  ind.missing <- which(Prob.Missing > max.miss.C & Prob.Missing <= max.missing.consider)
  out$InitialPvalues <- matrix(NA, nrow=p, ncol=K)
  Qvalues <- matrix(NA, nrow=p, ncol=K)
  max.H1 <- 0.5  #Maximum p-value to be considered H1. To be used if q-value fails
  for (g in ind.missing) {
    y.g <- Y[g,]
    out$InitialPvalues[g,] <- unlist(lapply(X = 1:K, function(k, y.g, Cov, C){ my.OLS(y = y.g[!is.na(y.g)], X = cbind(Cov,C[,k])[!is.na(y.g),])$p.value[ncol(Cov)+1] }, y.g=y.g, Cov=Cov, C=C))
  }
  for (k in 1:K) {
    ind.k <- !is.na(out$InitialPvalues[,k])
    p.k <- out$InitialPvalues[ind.k,k]
    try.k <- try(expr = {Q.k <- qvalue::qvalue(p.k); Qvalues[ind.k,k] <- Q.k$qvalues}, silent = TRUE)
    if (class(try.k) == "try-error") {
      n.0.k <- sum(p.k > max.H1)/(1-max.H1)
      Qvalues[ind.k,k] <- p.k*n.0.k/unlist(lapply(p.k, function(p.kg){sum(p.k<=p.kg)}))
    }
  }

  #  ###Preliminary estimate for missingness###
  #  n_cores <- max(detectCores() - 1, 1)
  #  cat(paste("Preliminary estimate of missingness mechanism using", n_cores, "cores..."))
  #  Exp.FI <- solve(t(out.C$W)%*%out.C$W)
  #  Tstat.start <- matrix(NA, p, K+1)
  #  
  #  cl <- makeCluster(n_cores)
  #  clusterExport(cl = cl, c("min.a", "max.a", "min.y0", "max.y0", "C", "W", "Exp.FI", "n.K.GMM", "t.df", "p.min.1", "p.min.2"), envir=environment())
  #  clusterEvalQ(cl = cl, expr = {source("/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/COPSAC_Metabolite/R/EstimateMissingnessGMM_t.R")})
  #  Tstat.start[ind.missing,] <- t( parApply(cl = cl, X = Y[ind.missing,], MARGIN = 1, function(y.g) {
  #    ind.g <- !is.na(y.g)
  #    miss.g <- Opt.GMM.t.grid(A.grid = seq(min.a,max.a,by=0.05), Pos.Grid = seq(min.y0,max.y0,by=0.2), y = y.g, C = C, K.ind = 1:n.K.GMM, n.iter = 2, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  #    prob.obs.g <- pt.gen(x = miss.g$a*y.g[ind.g]-miss.g$y0, df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
  #    prob.obs.g[prob.obs.g < 0.05] <- 0.05
  #    
  #    gamma.g <- solve(t(W[ind.g,]/prob.obs.g)%*%W[ind.g,], t(W[ind.g,]/prob.obs.g)%*%y.g[ind.g])
  #    resids.g <- W[ind.g,] * as.vector((y.g[ind.g] - W[ind.g,]%*%gamma.g))
  #    Var.g <- Exp.FI %*% (t(resids.g/prob.obs.g^2)%*%resids.g) %*% Exp.FI
  #    return(gamma.g / sqrt(diag(Var.g)))
  #  }) )
  #  stopCluster(cl)
  #
  #  cat("done\n")
  #  d <- ncol(out.C$W) - K
  #  Qvalues <- matrix(NA, nrow=p, ncol=K)
  #  out$InitialPvalues <- matrix(NA, nrow=p, ncol=K)
  #  max.H1 <- 0.5  #Maximum p-value to be considered H1. To be used if q-value fails
  #  for (k in 1:K) {
  #    ind.k <- !is.na(Tstat.start[,d+k])
  #    p.k <- 2*pnorm(-abs(Tstat.start[ind.k,d+k])); out$InitialPvalues[ind.k,k] <- p.k
  #    try.k <- try(expr = {Q.k <- qvalue::qvalue(p.k); Qvalues[ind.k,k] <- Q.k$qvalues}, silent = TRUE)
  #    if (class(try.k) == "try-error") {
  #      n.0.k <- sum(p.k > max.H1)/(1-max.H1)
  #      Qvalues[ind.k,k] <- p.k*n.0.k/unlist(lapply(p.k, function(p.kg){sum(p.k<=p.kg)}))
  #    }
  #  }
  
  
  ###Estimate missingness mechanism parameters with GMM###
  if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
  cat(paste("Estimating missingness mechanism using", n_cores, "cores..."))
  out$Theta.Miss <- matrix(NA, nrow=p, ncol=2)
  out$Theta.MAR <- matrix(NA, nrow=p, ncol=n.K.GMM+1)
  out$Pi.MAR <- matrix(NA, nrow=p, ncol=n)
  out$K.mat <- matrix(NA, nrow=p, ncol=n.K.GMM)
  out$Var.Theta <- vector(mode = "list", length = p)
  out$Q.use <- matrix(NA, nrow=p, ncol=n.K.GMM)
  out$Value.opt <- rep(NA, p)
  out$Pvalue.value <- rep(NA, p)
  out$Pvalue.value.gamma <- rep(NA, p)
  
  cl <- makeCluster(n_cores)
  clusterExport(cl = cl, c("Model.Pvalue", "min.a", "max.a", "min.y0", "max.y0", "C", "n.K.GMM", "t.df", "p.min.1", "p.min.2", "n.boot.J"), envir=environment())
  #clusterEvalQ(cl = cl, expr = {source("/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/COPSAC_Metabolite/R/EstimateMissingnessGMM_t.R")
  #                              source("/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/COPSAC_Metabolite/R/BayesianGMM.R")})
  clusterEvalQ(cl = cl, expr = {library(MetabMiss)})
  out.tmp <- parLapply( cl = cl, X = lapply(ind.missing, function(g){return(list(y=Y[g,], t.stat=Qvalues[g,], q=Qvalues[g,]))}), function(ll){
    out.par <- list()
    y.g <- ll$y
    #K.ind.g <- order(-abs(ll$t.stat))[1:n.K.GMM]
    K.ind.g <- order(ll$q)[1:n.K.GMM]
    
    #MAR fit#
    tmp <- glm.fit(x = Get.U(C = C, K.ind = K.ind.g), y = as.numeric(!is.na(ll$y)), family = binomial())
    out.par$theta.mar <- tmp$coefficients
    out.par$Pi.MAR <- tmp$fitted.values
    
    out.par$K.ind <- K.ind.g
    out.par$q <- ll$q[K.ind.g]
    
    out.2 <- Opt.GMM.t.grid(A.grid = seq(min.a,max.a,by=0.05), Pos.Grid = seq(min.y0,max.y0,by=0.2), y = y.g, C = C, K.ind = K.ind.g, n.iter = 2, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    out.par$value <- out.2$out.value*2/length(y.g)
    out.par$theta <- c(out.2$a, out.2$y0/out.2$a)
    out.par$value <- Val.Opt.t(a = out.par$theta[1], y0 = out.par$theta[2], y = y.g, C = C, K.ind = K.ind.g, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
    
    #Bootstrapped J test statistic#
    if (Model.Pvalue && n.K.GMM >= 2 && !(out.par$theta[1] <= (min.a+0.1) || out.par$theta[1] >= (max.a-0.1) || out.par$theta[2] <= (min.y0+0.2) || out.par$theta[2] >= (max.y0-0.2))) {
      tmp <- Boot.J(y = y.g, C = C, theta = out.par$theta, K.ind = K.ind.g, value.opt = out.par$value, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2, n.boot = n.boot.J)
      out.par$J.pvalue <- tmp$p.value
      if (!is.null(tmp$p.value.boot)) {
        out.par$J.pvalue.boot <- tmp$p.value.boot
      } else {
        out.par$J.pvalue.boot <- tmp$p.value
      }
    } else {
      out.par$J.pvalue <- NA
      out.par$J.pvalue.boot <- NA
    }
    
    #Variance#
    out.var.2 <- try(expr = {out.par$Var <- Var.GMM.t(theta = c(out.2$a, out.2$y0), y = y.g, C = C, K.ind = K.ind.g, W = out.2$W, t.df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)$Asy.Var}, silent = T)
    if (class(out.var.2) == "try-error") {
      out.par$Var <- NA
    }
    return(out.par)
  } )
  stopCluster(cl)
  out$Theta.Miss[ind.missing,] <- t(sapply(X = out.tmp, function(x){x$theta}))
  out$Theta.MAR[ind.missing,] <- t(sapply(X = out.tmp, function(x){x$theta.mar}))
  out$Pi.MAR[ind.missing,] <- t(sapply(X = out.tmp, function(x){x$Pi.MAR}))
  out$K.mat[ind.missing,] <- t(sapply(X = out.tmp, function(x){x$K.ind}))
  out$Var.Theta[ind.missing] <- lapply(X = out.tmp, function(x){x$Var})
  out$Q.use[ind.missing,] <- t(sapply(X = out.tmp, function(x){x$q}))
  out$Value.opt[ind.missing] <- as.vector(sapply(X = out.tmp, function(x){x$value}))
  out$Pvalue.value[ind.missing] <- unlist(lapply(X = out.tmp, function(x){x$J.pvalue.boot}))
  out$Pvalue.value.gamma[ind.missing] <- unlist(lapply(X = out.tmp, function(x){x$J.pvalue}))
  
  cat("done\n")
  if (n.K.GMM >= 2 && Model.Pvalue) {
    #tmp.q <- try(expr = {qvalue.model <- qvalue(out$Pvalue.value)})
    reject.BH <- BH.proc(p = out$Pvalue.value, alpha = BH.analyze.min)
    ind.q <- reject.BH == FALSE & !is.na(reject.BH)
    out$Ind.Confident <- ind.q
    out$BH.min <- BH.analyze.min
    #ind.q <- !is.na(out$Pvalue.value) & qvalue.model$lfdr >= lfdr.analyze.min
  } else {
    ind.q <- !is.na(out$Theta.Miss[,1]) & !(out$Theta.Miss[,1] <= (min.a+0.1) | out$Theta.Miss[,1] >= (max.a-0.1) | out$Theta.Miss[,2] <= (min.y0+0.2) | out$Theta.Miss[,2] >= (max.y0-0.2))
    out$Ind.Confident <- ind.q
  }
  
  
  ###Empirical Bayes to estimate prior for a and y0 and then MCMC to get posterior expectations and variances###
  if (shrink.Est) {
    Bayes.est <- match.arg(Bayes.est, choices = c("EmpBayes", "FullBayes", "FullBayes_ind"))
    out$Bayes.est <- Bayes.est
    
    #Empirical Bayes#
    Emp.Bayes.both <- Emp.Bayes.MuSigma.Both(Mu = cbind(log(out$Theta.Miss[ind.q,1]),out$Theta.Miss[ind.q,2]),
                                             Var = lapply( X = which(ind.q==T), function(g){tmp <- out$Var.Theta[[g]]; if (is.null(tmp)){(return(NA))}; diag(c(1/out$Theta.Miss[g,1],1))%*%tmp%*%diag(c(1/out$Theta.Miss[g,1],1))} ), simple.average=simple.average.EB)
    out$Emp.Bayes.loga <- Emp.Bayes.MuSigma(Mu.g = log(out$Theta.Miss[ind.q,1]), Var.g = unlist(lapply(X=out$Var.Theta[ind.q],function(x){x[1,1]}))/out$Theta.Miss[ind.q,1]^2, middle = "mean", shift.var = 0, refine.mu = T, simple.average=simple.average.EB)
    out$Emp.Bayes.y0 <- Emp.Bayes.MuSigma(Mu.g = out$Theta.Miss[ind.q,2], Var.g = unlist(lapply(X=out$Var.Theta[ind.q],function(x){x[2,2]})), middle = "mean", shift.var = 0, refine.mu = T, simple.average=simple.average.EB)
    out$V.prior <- Emp.Bayes.both$V
    out$mu.prior <- Emp.Bayes.both$mu
    
    out$Post.W <- matrix(1, nrow=p, ncol=n)
    out$Post.VarW <- matrix(1, nrow=p, ncol=n)
    out$Post.Theta <- matrix(NA, nrow=p, ncol=2)
    out$Post.Var <- vector(mode = "list", length = p)
    out$Post.Pi <- matrix(1, nrow=p, ncol=n)
    
    ind.q.MCMC <- !is.na(out$Theta.Miss[,1])
    
    if (Bayes.est == "FullBayes") {
      out.dep <- Bayes.GMM.Dependent(Y = Y[ind.q.MCMC,], C = C, K.ind = out$K.mat[ind.q.MCMC,], ind.variance = which(out$Ind.Confident[ind.q.MCMC] == T), p.min.1 = p.min.1, p.min.2 = p.min.2, t.df = t.df, y0.mean = Emp.Bayes.both$mu[2], log.a.mean = Emp.Bayes.both$mu[1], y0.start = out$Theta.Miss[ind.q.MCMC,2], a.start = out$Theta.Miss[ind.q.MCMC,1], n.iter = n.iter.MCMC, n.burn = n.burn.MCMC, prop.y0.sd = prop.y0.sd, prop.a.sd = prop.a.sd, include.norm = T, min.prob = min.prob.MCMC)
      out$Post.Theta[which(ind.q.MCMC==T),] <- out.dep$Post.Exp
      out$Post.Var[which(ind.q.MCMC==T)] <- out.dep$Post.Var
      out$Post.W[which(ind.q.MCMC==T),] <- out.dep$Post.Exp.W
      out$Post.VarW[which(ind.q.MCMC==T),] <- out.dep$Post.Var.W
      out$Post.Pi[which(ind.q.MCMC==T),] <- out.dep$Post.Exp.Pi
      out$Post.Exp.Var <- out.dep$Post.exp.Var
      out$Post.Var.Var <- out.dep$Post.Var.Var
    }
    
    if (Bayes.est == "FullBayes_ind") {
      out.indep <- Bayes.GMM.ind(Y = Y[ind.q.MCMC,], C = C, K.ind = out$K.mat[ind.q.MCMC,], ind.variance = which(out$Ind.Confident[ind.q.MCMC] == T), p.min.1 = p.min.1, p.min.2 = p.min.2, t.df = t.df, y0.mean = Emp.Bayes.both$mu[2], log.a.mean = Emp.Bayes.both$mu[1], y0.start = out$Theta.Miss[ind.q.MCMC,2], a.start = out$Theta.Miss[ind.q.MCMC,1], prop.y0.sd = prop.y0.sd, prop.a.sd = prop.a.sd, n.iter = n.iter.MCMC, n.burn = n.burn.MCMC, include.norm = T, min.prob = min.prob.MCMC)
      out$Post.Theta[which(ind.q.MCMC==T),] <- out.indep$Post.Exp
      out$Post.Var[which(ind.q.MCMC==T)] <- out.indep$Post.Var
      out$Post.W[which(ind.q.MCMC==T),] <- out.indep$Post.Exp.W
      out$Post.VarW[which(ind.q.MCMC==T),] <- out.indep$Post.Var.W
      out$Post.Pi[which(ind.q.MCMC==T),] <- out.indep$Post.Exp.Pi
      out$Post.Exp.Var <- out.indep$Post.exp.Var
      out$Post.Var.Var <- out$Post.Var.Var
    }
    
    if (Bayes.est == "EmpBayes") {
      #MCMC with covariance term = 0#
      if (is.null(n_cores)) {n_cores <- max(detectCores() - 1, 1)}
      cat(paste("MCMC using", n_cores, "cores..."))
      cl <- makeCluster(n_cores)
      clusterExport(cl = cl, c("C", "t.df", "p.min.1", "p.min.2", "Emp.Bayes.both", "prop.y0.sd", "prop.a.sd", "n.iter.MCMC", "n.burn.MCMC", "min.prob.MCMC"), envir=environment())
      #clusterEvalQ(cl = cl, expr = {source("/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/COPSAC_Metabolite/R/EstimateMissingnessGMM_t.R")
      #  source("/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/COPSAC_Metabolite/R/BayesianGMM.R")})
      clusterEvalQ(cl = cl, expr = {library(MetabMiss)})
      out.MCMC <- parLapply(cl = cl, X = lapply(which(ind.q.MCMC==T),function(g){return(list(K.ind=out$K.mat[g,], y=Y[g,], y0=out$Theta.Miss[g,2], a=out$Theta.Miss[g,1]))}), function(ll) {
        out.Bayes.g <- Bayes.GMM(y = ll$y, C = C, K.ind = ll$K.ind, p.min.1 = p.min.1, p.min.2 = p.min.2, t.df = t.df, y0.mean = Emp.Bayes.both$mu[2], log.a.mean = Emp.Bayes.both$mu[1], V = diag(diag(Emp.Bayes.both$V)), y0.start = ll$y0, a.start = ll$a, prop.y0.sd = prop.y0.sd, prop.a.sd = prop.a.sd, n.iter = n.iter.MCMC, n.burn = n.burn.MCMC, save.every = 0, include.norm = T, min.prob = min.prob.MCMC)
        return(list( theta=out.Bayes.g$Post.Exp, Var=out.Bayes.g$Post.Var, W=out.Bayes.g$Post.Exp.W, Pi=out.Bayes.g$Post.Exp.Pi, Var.W=out.Bayes.g$Post.Var.W ))
      })
      stopCluster(cl)
      cat("done\n")
      out$Post.Theta[which(ind.q.MCMC==T),] <- t(sapply(X = out.MCMC, function(x){x$theta}))
      out$Post.Var[which(ind.q.MCMC==T)] <- lapply(X = out.MCMC, function(x){x$Var})
      out$Post.W[which(ind.q.MCMC==T),] <- t(sapply(X = out.MCMC, function(x){x$W}))
      out$Post.VarW[which(ind.q.MCMC==T),] <- t(sapply(X = out.MCMC, function(x){x$Var.W}))
      out$Post.Pi[which(ind.q.MCMC==T),] <- t(sapply(X = out.MCMC, function(x){x$Pi}))
    }

  }
  out$ind.analyze <- sort(c(which(Prob.Missing <= max.miss.C), which(ind.q.MCMC == T)))
  return(out)
}


######Estimate coefficients of interest given missingness mechanism######
#This does not include confounder estimation

EstimateParam <- function(Y, X=NULL, Z=NULL, Miss.Obj, method=c("OLS", "Huber"), include.naive.est=F, Prob.Obs=NULL) {
  out <- list()
  method <- match.arg(method, c("OLS", "Huber"))
  Cov <- cbind(X,Z)
  
  if (is.null(X)) {
    d <- 0
  } else {
    X <- cbind(X)
    d <- ncol(X)
  }
  
  if (is.null(Z)) {
    r <- 0
  } else {
    Z <- cbind(Z)
    r <- ncol(Z)
  }
  
  n <- ncol(Y)
  p <- nrow(Y)
  ind.metab.missing <- which(!is.na(Miss.Obj$Theta.Miss[,1]) == T)
  ind.metab.obs <- Miss.Obj$ind.analyze[!(Miss.Obj$ind.analyze %in% ind.metab.missing)]
  
  #Get effect estimates#
  out$Beta <- matrix(NA, nrow=p, ncol=d+r)
  out$Beta.sw <- matrix(NA, nrow=p, ncol=d+r)   #Stabalized weights
  out$Beta.NoShrink <- matrix(NA, nrow=p, ncol=d+r)   #No shrinkage
  out$Var <- matrix(NA, nrow=p, ncol=d+r) 
  out$Var.sw <- matrix(NA, nrow=p, ncol=d+r); out$Var.sw.vc <- matrix(NA, nrow=p, ncol=d+r)
  out$Var.NoShrink <- matrix(NA, nrow=p, ncol=d+r) 
  if (include.naive.est) {
    out$Beta.naive <- matrix(NA, nrow=p, ncol=d+r)
    out$Var.naive <- matrix(NA, nrow=p, ncol=d+r)
    tmp <- t(apply(X = Y[Miss.Obj$ind.analyze,], MARGIN = 1, function(y) {
      ind.obs <- !is.na(y)
      beta <- solve(t(Cov[ind.obs,])%*%Cov[ind.obs,],t(Cov[ind.obs,])%*%y[ind.obs])
      sigma2 <- sum((y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))^2)/(sum(ind.obs)-d-r)
      return(c(beta,sigma2*diag(solve(t(Cov[ind.obs,])%*%Cov[ind.obs,]))))
    }))
    out$Beta.naive[Miss.Obj$ind.analyze,] <- tmp[,1:(d+r)]
    out$Var.naive[Miss.Obj$ind.analyze,] <- tmp[,(d+r+1):(2*d+2*r)]
  }
  
  if (method == "OLS") {
    #Fully observed metabolites#
    tmp <- t(apply(X = Y[ind.metab.obs,], MARGIN = 1, function(y) {
      ind.obs <- !is.na(y)
      beta <- solve(t(Cov[ind.obs,])%*%Cov[ind.obs,],t(Cov[ind.obs,])%*%y[ind.obs])
      sigma2 <- sum((y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))^2)/(n-d-r)
      return(c(beta,sigma2*diag(solve(t(Cov[ind.obs,])%*%Cov[ind.obs,]))))
    }))
    out$Beta[ind.metab.obs,] <- tmp[,1:(d+r)]
    out$Beta.sw[ind.metab.obs,] <- tmp[,1:(d+r)]
    out$Beta.NoShrink[ind.metab.obs,] <- tmp[,1:(d+r)]
    out$Var[ind.metab.obs,] <- tmp[,(d+r+1):(2*d+2*r)]
    out$Var.sw[ind.metab.obs,] <- tmp[,(d+r+1):(2*d+2*r)]; out$Var.sw.vc[ind.metab.obs,] <- tmp[,(d+r+1):(2*d+2*r)]
    out$Var.NoShrink[ind.metab.obs,] <- tmp[,(d+r+1):(2*d+2*r)]
    
    #Partially observed metabolites#
    tmp <- t(sapply(X = ind.metab.missing, function(g) {
      y <- Y[g,]
      ind.obs <- !is.na(y)
      weights <- Miss.Obj$Post.W[g,ind.obs]; var.weights <- Miss.Obj$Post.VarW[g,ind.obs]
      beta <- solve( t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,], t(Cov[ind.obs,]*weights)%*%y[ind.obs] )
      score <- Cov[ind.obs,]*(y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))
      Sample.FI <- solve(t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,])
      return( c(beta,diag(Sample.FI%*%( t(score*(weights^2+var.weights))%*%score )%*%Sample.FI)) )
    }))
    out$Beta[ind.metab.missing,] <- tmp[,1:(d+r)]
    out$Var[ind.metab.missing,] <- tmp[,(d+r+1):(2*d+2*r)]
    
    tmp <- t(sapply(X = ind.metab.missing, function(g) {
      y <- Y[g,]
      ind.obs <- !is.na(y)
      weights <- Miss.Obj$Post.W[g,ind.obs]*Miss.Obj$Pi.MAR[g,ind.obs]; var.weights <- Miss.Obj$Post.VarW[g,ind.obs]*Miss.Obj$Pi.MAR[g,ind.obs]^2
      beta <- solve( t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,], t(Cov[ind.obs,]*weights)%*%y[ind.obs] )
      score <- Cov[ind.obs,]*(y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))  #n x (d+r)
      Sample.FI <- solve(t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,])
      
      cov.tmp <- cbind(Cov[ind.obs,])
      #score.vc <- t(sapply(X = 1:nrow(score), function(i){ as.vector(solve( diag(1,ncol(cov.tmp),ncol(cov.tmp)) - weights[i]*cbind(cov.tmp[i,])%*%rbind(cov.tmp[i,])%*%Sample.FI, score[i,] )) }))
      Weights.vc <- 1/(1 - weights*apply(cov.tmp, 1, function(x) { sum(x*as.vector(Sample.FI%*%x)) }))
      
      return( c(beta,diag(Sample.FI%*%( t(score*(weights^2+var.weights))%*%score )%*%Sample.FI),diag(Sample.FI%*%( t(score*Weights.vc^2*(weights^2+var.weights))%*%score )%*%Sample.FI)) )
    }))
    out$Beta.sw[ind.metab.missing,] <- tmp[,1:(d+r)]
    out$Var.sw[ind.metab.missing,] <- tmp[,(d+r+1):(2*d+2*r)]; out$Var.sw.vc[ind.metab.missing,] <- tmp[,(2*d+2*r+1):(3*d+3*r)]
    
    tmp <- t(sapply(X = ind.metab.missing, function(g) {
      y <- Y[g,]
      ind.obs <- !is.na(y)
      weights <- 1/pt.gen(x = Miss.Obj$Theta.Miss[g,1]*(y[ind.obs]-Miss.Obj$Theta.Miss[g,2]), df = t.df, p.min.1 = p.min.1, p.min.2 = p.min.2)
      beta <- solve( t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,], t(Cov[ind.obs,]*weights)%*%y[ind.obs] )
      score <- Cov[ind.obs,]*(y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))
      Sample.FI <- solve(t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,])
      return( c(beta,diag(Sample.FI%*%( t(score*weights^2)%*%score )%*%Sample.FI)) )
    }))
    out$Beta.NoShrink[ind.metab.missing,] <- tmp[,1:(d+r)]
    out$Var.NoShrink[ind.metab.missing,] <- tmp[,(d+r+1):(2*d+2*r)]
    
    if (!is.null(Prob.Obs)) {
      out$Beta.obs <- matrix(NA, nrow=p, ncol=d+r)
      out$Var.obs <- matrix(NA, nrow=p, ncol=d+r)
      tmp <- t(sapply(X = ind.metab.missing, function(g) {
        y <- Y[g,]
        ind.obs <- !is.na(y)
        weights <- 1/Prob.Obs[g,ind.obs]
        beta <- solve( t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,], t(Cov[ind.obs,]*weights)%*%y[ind.obs] )
        score <- Cov[ind.obs,]*(y[ind.obs]-as.vector(Cov[ind.obs,]%*%beta))
        Sample.FI <- solve(t(Cov[ind.obs,]*weights)%*%Cov[ind.obs,])
        return( c(beta,diag(Sample.FI%*%( t(score*weights^2)%*%score )%*%Sample.FI)) )
      }))
      out$Beta.obs[ind.metab.missing,] <- tmp[,1:(d+r)]
      out$Var.obs[ind.metab.missing,] <- tmp[,(d+r+1):(2*d+2*r)]
    }
  }
  
  return(out)
}


####OLS with weights####
Standard.OLS <- function(y, Cov, weights=NULL, Cov.total=NULL) {
  Cov <- cbind(Cov)
  n <- length(y)
  
  if (is.null(weights)) {
    b <- solve(t(Cov)%*%Cov,t(Cov)%*%y)
    sigma2 <- 1/(n-ncol(Cov))*sum((y-as.vector(Cov%*%b))^2)
    return(list(b=b, z=b/sqrt(sigma2)/sqrt(diag(solve(t(Cov)%*%Cov)))))
  }
  b <- as.vector(solve(t(Cov*weights)%*%Cov,t(Cov*weights)%*%y))
  hess <- solve(t(Cov*weights)%*%Cov)
  Inflate <- 1/(1 - weights * apply(X = Cov, MARGIN = 1, function(x){sum(x*as.vector(hess%*%x))}))
  Mid <- t(Cov*weights^2*Inflate^2*(y - as.vector(Cov%*%b))^2)%*%Cov
  if (is.null(Cov.total)) {
    return(list(b=b, z=b/sqrt(diag(hess%*%Mid%*%hess))))
  }
  hess <- solve(t(Cov.total)%*%Cov.total)
  return(list(b=b, z=b/sqrt(diag(hess%*%Mid%*%hess))))
} 


#####BH-FDR#####

BH.proc <- function(p, alpha=0.25) {
  n <- length(p)
  Reject <- rep(NA, n)   #Reject null?
  p.obs <- p[!is.na(p)]
  n.obs <- length(p.obs)
  
  Reject.obs <- rep(FALSE, n.obs)
  order.obs <- order(p.obs)
  ind.less <- which(p.obs[order.obs] - (1:n.obs)/n.obs*alpha < 0)
  if (length(ind.less) > 0) {
    Reject.obs[p.obs <= max(p.obs[order.obs][ind.less])] <- TRUE
  }
  Reject[!is.na(p)] <- Reject.obs
  return(Reject)
}

######OLS######

my.OLS <- function(y, X, d.add = 0) {
  n <- length(y)
  X <- cbind(X)
  d <- ncol(X)
  hess <- solve(t(X)%*%X)
  beta.hat <- hess%*%t(X)%*%y
  resids <- as.vector(y - X %*% beta.hat)
  sigma2.hat <- sum( resids^2 ) / (n - d - d.add)
  var.inflated.beta <- hess%*%(t(X * (resids^2 / ( 1-rowSums((X%*%hess)*X) )^2)) %*% X)%*%hess
  return(list(beta=beta.hat, var.beta.hat=sigma2.hat*diag(hess), var.inflated.beta=diag(var.inflated.beta), sigma2=sigma2.hat, p.value=2*pnorm(-abs( beta.hat / sqrt(sigma2.hat) / sqrt(diag( hess )) ))))
}

######Choose the number of potential instruments######

#' Choose the number of potential instruments. This is \code{K} in the function \code{EstimateMissing}.
#'
#' Choose the number of potential instruments using a q-value threshold
#' 
#' @param Y a \code{p} x \code{n} data matrix of log2-transformed metabolite intensities, where \code{p} = #of metabolites and \code{n} = #of samples. Missing values should be left as \code{NA}. This is the only variable that must be specified.
#' @param Cov a \code{n} x \code{d} matrix of covariates, where d <= 2. The default, and recommended value, is the vector of all 1s.
#' @param max.miss.C a number between 0 and 1. The maximum fraction of missing data a metabolite is allowed to have to be considered nearly completely observed. The default, and recommended value, is 0.05.
#' @param max.missing.consider a number between 0 and 1. The maximum fraction of missing data a metabolite is allowed to have. Any metabolites with missingness fractions greater than this will be ignored. Default, and recommended value, is 0.5.
#' @param K.max an integer >=2. The maximum number of potential instruments to consider. If unspecified, it is set to be sva::num.sv estimate for K applied to the metabolites with complete data.
#' @param q.thresh a vector of numbers between 0 and 1. The q-value thresholds to consider. It defaults to c(0.01, 0.05, 0.1).
#' 
#' @return a list. \item{Frac1}{a \code{K.max} x \code{length(q.thresh)} matrix. The (k,j)th entry is the fraction of metabolites with at least one factor 1,...,k with q-value less than or equal to \code{q.thresh[j]}.} \item{Frac2}{a \code{K.max} x \code{length(q.thresh)} matrix. The (k,j)th entry is the fraction of metabolites with at least two factors 1,...,k with q-value less than or equal to \code{q.thresh[j]}. It is returned only if \code{d = 1}} \item{frac.thresh}{The fraction of factors with q-value < 0.05. If this is < 0.75, we recommend not estimating the missingness mechanism.} \item{K.recommended}{The recommended K (a scalar). This can be used as \code{K} in the function \code{EstimateMissing}. This is the smallest number of factors such that \code{frac.thresh} of the metabolites with missing data have at least 1 factor with q-value < 0.05.}
#' 
#' @export
Num.Instruments <- function(Y, Cov=NULL, max.miss.C = 0.05, max.missing.consider=0.5, K.max=NULL, q.thresh=c(0.01, 0.05, 0.1)) {
  q.recommended <- 0.05
  min.thresh <- 0.9
  p <- nrow(Y)
  n <- ncol(Y)
  if (is.null(Cov)) {Cov <- rep(1,n)}
  Cov <- cbind(Cov); d <- ncol(Cov)
  Frac.Missing <- rowMeans(is.na(Y))
  if (is.null(K.max)) {K.max <- sva::num.sv(dat = Y[Frac.Missing==0,], mod = Cov)}
  N1.thresh <- matrix(NA, nrow=K.max, ncol=length(q.thresh))
  N1.recommended <- rep(NA, K.max)
  if (d == 1) {N2.thresh <- matrix(NA, nrow=K.max, ncol=length(q.thresh))}
  
  max.H1 <- 0.5  #Maximum p-value to be considered H1. To be used if q-value fails
  ind.missing <- which((Frac.Missing > max.miss.C & Frac.Missing <= max.missing.consider) == T)
  for (K in 2:K.max) {
    cat(paste0("K = ", K, "/", K.max, "..."))
    C <- EstC.0(Y = Y, K = K, Cov = Cov, max.miss = max.miss.C, max.iter = 1000, n.repeat.Sigma = 1)$C
    
    Q.K <- matrix(NA, nrow=p, ncol=K)
    InitialPvalues <- matrix(NA, nrow=p, ncol=K)
    for (g in ind.missing) {
      y.g <- Y[g,]
      InitialPvalues[g,] <- unlist(lapply(X = 1:K, function(k, y.g, Cov, C){ my.OLS(y = y.g[!is.na(y.g)], X = cbind(Cov,C[,k])[!is.na(y.g),])$p.value[ncol(Cov)+1] }, y.g=y.g, Cov=Cov, C=C))
    }
    for (k in 1:K) {
      ind.k <- !is.na(InitialPvalues[,k])
      p.k <- InitialPvalues[ind.k,k]
      try.k <- try(expr = {Q.k <- qvalue::qvalue(p.k); Q.K[ind.k,k] <- Q.k$qvalues}, silent = TRUE)
      if (class(try.k) == "try-error") {
        n.0.k <- sum(p.k > max.H1)/(1-max.H1)
        Q.K[ind.k,k] <- p.k*n.0.k/unlist(lapply(p.k, function(p.kg){sum(p.k<=p.kg)}))
      }
    }
    
    Min.Q.K <- apply(X = Q.K, MARGIN = 1, min)
    for (j in 1:length(q.thresh)) {
      N1.thresh[K,j] <- sum(Min.Q.K <= q.thresh[j], na.rm = T) / length(ind.missing)
    }
    N1.recommended[K] <- sum(Min.Q.K <= q.recommended, na.rm = T) / length(ind.missing)
    if (d == 1) {
      Min2.Q.K <- apply(X = Q.K, MARGIN = 1, function(x){sort(x)[2]})
      for (j in 1:length(q.thresh)) {
        N2.thresh[K,j] <- sum(Min2.Q.K <= q.thresh[j], na.rm = T) / length(ind.missing)
      }
    }
    cat("done\n")
  }
  tmp.recommend <- which(N1.recommended >= min.thresh)
  while (length(tmp.recommend) == 0 && min.thresh > 0) {
    min.thresh <- max(0, min.thresh - 0.05)
    tmp.recommend <- which(N1.recommended >= min.thresh)
  }
  K.recommended <- max(2, min(tmp.recommend))
  
  rownames(N1.thresh) <- 1:K.max
  colnames(N1.thresh) <- q.thresh
  if (d == 1) {
    rownames(N2.thresh) <- 1:K.max
    colnames(N2.thresh) <- q.thresh
    return(list(Frac1=N1.thresh, Frac2=N2.thresh, K.recommended=K.recommended, frac.thresh=min.thresh))
  }
  return(list(Frac1=N1.thresh, K.recommended=K.recommended, frac.thresh=min.thresh))
}

