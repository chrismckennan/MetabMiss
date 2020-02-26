##Modified quantile normalization##
my.quant <- function(X, min.quant=NULL) { #p x n data matrix
  p <- nrow(X)
  n <- ncol(X)
  p.missing <- apply(X, 2, function(x){sum(is.na(x))})   #Number of missing analytes in each sample (an n-vector)
  if (is.null(min.quant)) {min.quant <- max(p.missing) + 1}    #Everything below this is standardized in the same way
  Quants <- apply(X, 2, function(x){ x[is.na(x)] <- 0; quantile(x, (min.quant:p)/p) })
  mean.quants <- apply(Quants , 1, mean)  #mean of quantiles. The first element is the in.quant/p quantile
  
  out <- X
  for (i in 1:n) {
    q.i <- Quants[,i]
    x.i <- out[,i]
    x.order <- x.i; x.order[is.na(x.order)] <- 0; order.i <- rank(x.order, ties.method = "first")
    out[order.i < min.quant,i] <- x.i[order.i < min.quant] + mean.quants[1] - q.i[1]
    
    x.rest <- x.i[order.i >= min.quant]
    order.rest <- rank(x.rest, ties.method = "first")
    out[order.i >= min.quant,i] <- mean.quants[order.rest]
  }
  return(out)
}


####Enrichment p-values####
#Performs a simple chi-squared test for independence

GetTable <- function(x.rows, x.cols, t.rows, t.cols, t.rows.type="greater", t.cols.type) {
  ind.remove <- is.na(x.rows) | is.na(x.cols)
  x.rows <- x.rows[!ind.remove]; x.cols <- x.cols[!ind.remove]
  if (t.rows.type == "greater") {
    ind.rows <- x.rows >= t.rows
  } else {
    ind.rows <- x.ros <= t.rows
  }
  
  if (t.cols.type == "greater") {
    ind.cols <- x.cols >= t.cols
  } else {
    ind.cols <- x.cols <= t.cols
  }
  
  return( Enrichment(n.ss=sum(ind.rows & ind.cols), n.sn=sum(ind.rows & !ind.cols), n.ns=sum(!ind.rows & ind.cols), n.nn=sum(!ind.rows & !ind.cols)) )
}

Enrichment.indices <- function(ind.rows, ind.cols, log.p = T) {
  ind.non.na <- !(is.na(ind.rows) | is.na(ind.cols))
  Enrichment(n.ss=sum(ind.rows & ind.cols & ind.non.na), n.sn=sum(ind.rows & !ind.cols & ind.non.na), n.ns=sum(!ind.rows & ind.cols & ind.non.na), n.nn=sum(!ind.rows & !ind.cols & ind.non.na), log.p = log.p)
}

Enrichment <- function(n.ss, n.sn, n.ns, n.nn, log.p = T) {
  n <- n.ss + n.sn + n.ns + n.nn
  p.ss <- n.ss/n
  p.ns <- n.ns/n
  p.sn <- n.sn/n
  p.nn <- n.nn/n
  
  p.s0 <- (n.ss + n.sn)/n; p.n0 <- 1-p.s0
  p.0s <- (n.ss + n.ns)/n; p.0n <- 1-p.0s
  t <- n * ((p.ss-p.s0*p.0s)^2/(p.s0*p.0s) + (p.ns-p.n0*p.0s)^2/(p.n0*p.0s) + (p.sn-p.s0*p.0n)^2/(p.s0*p.0n) + (p.nn-p.n0*p.0n)^2/(p.n0*p.0n))
  return(pchisq(t, df=1, lower.tail = F, log.p = log.p))
}

###Ordinary least squares###

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

###Multivariate linear model###
#X is assumed to be univariate
Multivariate.OLS <- function(Y, X, Z, missing.max=0.5) {
  Prob.Miss <- apply(Y, 1, function(y){mean(is.na(y))})
  p <- nrow(Y)
  X <- cbind(X)
  Cov <- cbind(X,Z)
  
  Beta <- rep(NA, p)
  Var.Beta <- rep(NA, p)
  Var.Beta.inflated <- rep(NA, p)
  Sigma <- rep(NA, p)
  p.Beta <- rep(NA, p)
  
  Beta.all <- Y[Prob.Miss==0,]%*%Cov%*%solve(t(Cov)%*%Cov)
  Beta[Prob.Miss==0] <- Beta.all[,1]
  Sigma[Prob.Miss==0] <- rowSums((Y[Prob.Miss==0,]-Beta.all%*%t(Cov))^2)/(nrow(X)-ncol(Cov))
  Var.Beta[Prob.Miss==0] <- solve(t(Cov)%*%Cov)[1,1] * Sigma[Prob.Miss==0]
  Var.Beta.inflated[Prob.Miss==0] <- Var.Beta[Prob.Miss==0]
  p.Beta[Prob.Miss==0] <- 2*pnorm(-abs(Beta[Prob.Miss==0]) / sqrt(Var.Beta[Prob.Miss==0]))
  
  ind.miss <- which(Prob.Miss > 0 & Prob.Miss <= missing.max)
  tmp <- lapply(ind.miss, function(g){ ind.g <- !is.na(Y[g,]); my.OLS(y = Y[g,ind.g], X = Cov[ind.g,]) })
  Beta[ind.miss] <- unlist(lapply(tmp,function(x){x$beta[1]}))
  Sigma[ind.miss] <- unlist(lapply(tmp,function(x){x$sigma2}))
  Var.Beta[ind.miss] <- unlist(lapply(tmp,function(x){x$var.beta.hat[1]}))
  Var.Beta.inflated[ind.miss] <- unlist(lapply(tmp,function(x){x$var.inflated.beta[1]}))
  p.Beta[ind.miss] <- unlist(lapply(tmp,function(x){x$p.value[1]}))
  
  return(list(B=Beta, Sigma=Sigma, Var.Beta=Var.Beta, Var.Beta.inflated=Var.Beta.inflated, p.Beta=p.Beta))
}

###permutation pvalue###

my.perm.pvalue <- function(all.data, data, n.boot = 500) {
  n <- length(data)/2
  expected <- get.expected(all.data, n=n)
  actual <- get.actual(data)
  stat.chisq <- sum((actual - expected)^2/expected)
  p.both <- sum(expected[1:2])/n
  stat.both <- (sum(actual[1:2])/n-p.both)/sqrt(p.both*(1-p.both)/n)
  
  #stat.boot <- replicate(n.boot, {data.boot <- sample(all.data, size = 2*n, replace=F); expected.boot <- get.expected(data.boot); sum((get.actual(data.boot) - expected.boot)^2/expected.boot)})
  #return( c((sum(stat.boot >= stat) + 1)/(n.boot + 1), 1-pchisq(stat, df=2)) )
  c(1-pchisq(stat.chisq, df=2), 1-pnorm(stat.both))
}

get.actual <- function(data) {
  data.na <- matrix(is.na(data), nrow = length(data)/2, ncol=2, byrow = TRUE)
  c(sum(apply(data.na, 1, prod)), sum(apply(!data.na, 1, prod)), sum(apply(data.na, 1, function(x){ (x[1] & !x[2]) | (!x[1] & x[2]) })))
}

get.expected <- function(data, n) {
  p.g <- sum(is.na(data))/length(data)
  c(n*p.g^2, n*(1-p.g)^2, 2*n*p.g*(1-p.g))
}

####Compute retention time break points####
breaks.rt <- function(bin.num, rt, delta.ends = 0.01) {
  ids <- sort(unique(bin.num))
  break.points <- rep(NA, length(ids)+1)
  for (i in 1:length(ids)) {
    if (i == 1) {
      break.points[1:2] <- c(min(rt[bin.num==ids[i]]) - delta.ends, min(rt[bin.num==ids[i+1]]) - delta.ends)
    } else {
      if (i == length(ids)) {
        break.points[i+1] <- max(rt) + delta.ends
      } else {
        break.points[i+1] <- min(rt[bin.num==ids[i+1]]) - delta.ends
      }
    }
  }
  return(break.points)
}

####which function for indices of a matrix####

Which.max <- function(X, MARGIN=1) {
  p <- nrow(X)
  n <- ncol(X)
  max.margin <- apply(X, MARGIN = MARGIN, function(x){max(x, na.rm = T)})
  ind.margin <- which.max(max.margin)
  if (MARGIN==1) {
    return( c( ind.margin, which.max(X[ind.margin,], na.rm=T) ) )
  } else {
    return( c( which.max(X[,ind.margin], na.rm=T), ind.margin ) )
  }
}

Which.min <- function(X, MARGIN=1) {
  Which.max(-X, MARGIN=MARGIN)
}

###Orthogonal projections###
Compute.Q <- function(X) {
  X <- cbind(X)
  t <- qr(X)
  return( qr.Q(t, complete = T)[,(t$rank+1):nrow(X)] )
}

###Skew###
skew <- function(x) {
  x <- x[!is.na(x)]
  x <- scale(x, center = T, scale = T)
  mean(x^3)
}


###Pvalue plot###

Pvalue.plot <- function(pvalues, use.log=T, xlab=NULL, ylab=NULL, col.points=0, col.indices=NULL, col="blue", pch.col=1, lwd.abline=1, cex.indices=1, ...) {
  n <- sum(!is.na(pvalues))
  if (!is.null(col.indices)) {col.indices <- match(x = col.indices, table = which(!is.na(pvalues)==T))}
  pvalues <- pvalues[!is.na(pvalues)]
  seq.plot <- (1:n)/(n+1)
  if (use.log) {
    if (is.null(xlab)) {xlab <- "-log10(Expected P-value)"}
    if (is.null(ylab)) {ylab <- "-log10(P-value)"}
    if (!is.null(col.indices)) {
      plot(-log10(seq.plot), -log10(sort(pvalues)), pch="", xlab=xlab, ylab=ylab, ...)
      points( -log10(rank(pvalues)[-col.indices]/(n+1)), -log10(pvalues[-col.indices]) )
      points( -log10(rank(pvalues)[col.indices]/(n+1)), -log10(pvalues[col.indices]), col=col, pch=pch.col, cex=cex.indices)
    } else {
      if (col.points > 0) {
        plot(-log10(seq.plot), -log10(sort(pvalues)), pch="", xlab=xlab, ylab=ylab, ...)
        points( -log10(seq.plot[(col.points+1):n]), -log10(sort(pvalues)[(col.points+1):n]) )
        points( -log10(seq.plot[1:col.points]), -log10(sort(pvalues)[1:col.points]), col=col, pch=pch.col)
      } else {
        plot(-log10(seq.plot), -log10(sort(pvalues)), xlab=xlab, ylab=ylab, ...)
      }
    }
    abline(a=0,b=1,col="red",lwd=lwd.abline)
  } else {
    if (is.null(xlab)) {xlab <- "Expected P-value"}
    if (is.null(ylab)) {ylab <- "P-value"}
    plot(seq.plot, sort(pvalues), xlab=xlab, ylab=ylab, ...)
    abline(a=0,b=1,col="red",lwd=lwd.abline)
  }
}


###Big boxplot###
Big.boxplot <- function(sim.results, q.plot=seq(0.05, 0.25, by=0.05), i.max=NULL, names.plot=NULL, width.split=1) {
  names.results <- names(sim.results)
  if (is.null(names.plot)) {names.plot <- names.results}
  if (is.null(i.max)) {i.max <- length(sim.results[[ names.results[1] ]][,1])}
  
  for (i in 1:length(q.plot)) {
    for (j in 1:length(names.results)) {
      if (i == 1 && j == 1) {
        mat <- cbind( sim.results[[ names.results[[j]] ]][1:i.max,i] )
      } else {
        mat <- cbind(mat, sim.results[[ names.results[[j]] ]][1:i.max,i])
      }
      if (j == length(names.results) && i < length(q.plot)) {
        mat <- cbind(mat, matrix(-1, nrow=i.max, ncol=width.split))
      }
    }
  }
  
  n.q <- length(q.plot)
  n.names <- length(names.results)
  boxplot(mat, axes=F, at=1:ncol(mat), use.cols = T, ylim=c(0,max(mat)))
  axis(side = 2, at = seq(0,1,by=0.05))
  axis(side = 1, at = 1:ncol(mat), labels = rep(c(names.plot, rep("",length=width.split)), times=length(q.plot))[-((ncol(mat)+1):(n.q*(n.names+width.split)))], las = 2)
  for (i in 1:length(q.plot)) {
    min.point <- (i-1)*(n.names+width.split) + 1
    max.point <- (i-1)*(n.names+width.split) + n.names
    lines(seq(min.point, max.point, length=100), rep(q.plot[i], length=100), col="red", lty=1, lwd=2)
  }
}

##Genome inflation factor##
Inflation.factor <- function(z) {
  z <- z[!is.na(z)]
  if (sum(z >= 0 & z <= 1) == length(z)) {return( c(median(z)/0.5, mean(z)/0.5) )}
  c(median(abs(z)) / sqrt(2/pi), sqrt(mean(z^2)))
}

###Confidence intervals###
Compute.CI <- function(Z, B, P.miss, max.i, min.B, max.B, min.miss=0.05, max.miss=0.5, alpha=0.05) {
  tmp <- unlist(lapply(1:max.i,function(i){ind.i<-abs(B[,i])>min.B&abs(B[,i])<=max.B&P.miss[,i]>min.miss&P.miss[,i]<=max.miss; Z[ind.i,i]>=qnorm(alpha/2)&Z[ind.i,i]<=qnorm(1-alpha/2)}))
  return(list(cov=mean(tmp), n=length(tmp)))
}

###GIF###
GIF.lambda <- function(z) {
  z <- z[!is.na(z)]
  if (sum(z <= 1 & z >= 0) == length(z)) {
    return(qnorm(median(z)/2)^2/0.456)
  }
  return(median(z^2)/0.456)
}

Get.RMSE <- function(beta, beta.hat, ind.use) {
  beta <- beta[ind.use]; beta.hat <- beta.hat[ind.use]
  sqrt(mean((beta[!is.na(beta.hat)]-beta.hat[!is.na(beta.hat)])^2))
}

###Q###
Get.Q <- function(Y, C, Cov=NULL, min.miss=0.05, max.miss=0.5) {
  max.H1 <- 0.5
  C <- cbind(C)
  n <- nrow(C)
  if (is.null(Cov)) {Cov <- rep(1,n)}
  Cov <- cbind(Cov)
  K <- ncol(C)
  frac.miss <- rowMeans(is.na(Y))
  ind.missing <- which(frac.miss > min.miss & frac.miss <= max.miss)
  Q.K <- matrix(NA, nrow=p, ncol=K)
  InitialPvalues <- matrix(NA, nrow = p, ncol = K)
  for (g in ind.missing) {
    y.g <- Y[g, ]
    InitialPvalues[g, ] <- unlist(lapply(X = 1:K, function(k, y.g, Cov, C) {
                my.OLS(y = y.g[!is.na(y.g)], X = cbind(Cov, C[,k])[!is.na(y.g), ])$p.value[ncol(Cov) + 1]
                 }, y.g = y.g, Cov = Cov, C = C))
  }
  for (k in 1:K) {
    ind.k <- !is.na(InitialPvalues[, k])
    p.k <- InitialPvalues[ind.k, k]
    try.k <- try(expr = {
      Q.k <- qvalue::qvalue(p.k)
      Q.K[ind.k, k] <- Q.k$qvalues
    }, silent = TRUE)
    if (class(try.k) == "try-error") {
      n.0.k <- sum(p.k > max.H1)/(1 - max.H1)
      Q.K[ind.k, k] <- p.k * n.0.k/unlist(lapply(p.k, function(p.kg){sum(p.k <= p.kg)}))
    }
  }
  return(list(Q=Q.K,p=InitialPvalues))
}

