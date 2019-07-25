
###################################################### 
# High-Dimensional Cumulative Covariance Estimator 
# by Factor Modeling and Regularization
######################################################

cce.factor <- function(yuima, method = "HY", factor = NULL, PCA = FALSE, 
                       nfactor = "interactive", regularize = "glasso",
                       taper, group = 1:(dim(yuima) - length(factor)),
                       lambda = "bic", weight = TRUE, nlambda = 10, 
                       ratio, N, thr.type = "soft", thr = NULL, 
                       tau = NULL, par.alasso = 1, par.scad = 3.7, 
                       thr.delta = 0.01, frequency = 300, utime, ...){
  
  cmat <- cce(yuima, method = method, frequency = frequency,utime = utime, ...)$covmat
  
  if(missing(N)) N <- effective.n(yuima, method, frequency, utime)
  
  if(PCA){
    
    ed <- eigen(cmat, symmetric = TRUE)
    pc <- ed$values
    
    if(nfactor == "interactive"){
      plot(pc, type = "h", ylab = "PC scores", main = "Scree plot")
      nfactor <- readline("Type the number of factors you want to use. \n")
    }
    
    if(nfactor == 0){
      
      factor <- NULL
      
      sigma.z <- sigma.y
      beta.hat <- NULL
      sigma.x <- NULL
      
    }else{
      
      factor <- "PCA"
      
      sigma.x <- diag(ed$values[1:nfactor], nfactor)
      inv.sigma.x <- solve(sigma.x)
      beta.hat <- as.matrix(ed$vectors[ ,1:nfactor])
      sigma.f <- tcrossprod(beta.hat %*% sigma.x, beta.hat)
      sigma.z <- cmat - sigma.f
      
    }
    
  }else{
    
    pc <- NULL
    
    if(is.null(factor)){
      
      sigma.z <- cmat
      beta.hat <- NULL
      sigma.x <- NULL
      
    }else{
      
      if(is.character(factor))
        factor <- which(rownames(cmat) %in% factor)
      
      sigma.y <- as.matrix(cmat[-factor,-factor])
      sigma.x <- as.matrix(cmat[factor,factor])
      sigma.yx <- as.matrix(cmat[-factor,factor])
      
      # is sigma.x singular?
      
      inv.sigma.x <- solve(sigma.x)
      beta.hat <- sigma.yx %*% inv.sigma.x
      sigma.f <- tcrossprod(beta.hat %*% sigma.x, beta.hat)
      sigma.z <- sigma.y - sigma.f
      
    }
    
  }
  
  reg.result <- cce.regularize(sigma.z, regularize, taper, group,
                               lambda, weight, nlambda, ratio, N,
                               thr.type, thr, tau, par.alasso, par.scad, thr.delta)
  
  if(is.null(factor)){
    covmat.y <- reg.result$cmat
    premat.y <- reg.result$pmat
  }else{
    covmat.y <- sigma.f + reg.result$cmat
    # using Sherman-Morisson-Woodbury formula to compute the inverse
    beta.pmat <- crossprod(beta.hat, reg.result$pmat)
    premat.y <- reg.result$pmat - 
      crossprod(beta.pmat, 
                solve(inv.sigma.x + beta.pmat %*% beta.hat) %*% beta.pmat)
  }
  
  return(list(covmat.y = covmat.y, premat.y = premat.y,
              beta.hat = beta.hat, covmat.x = sigma.x,
              covmat.z = reg.result$cmat, 
              premat.z = reg.result$pmat, 
              sigma.z = sigma.z, pc = pc))
}


cce.regularize <- function(cmat, method = "tapering", taper, group,
                           lambda = "bic", weight = TRUE, nlambda = 10, 
                           ratio, N, thr.type = "hard", thr = NULL, 
                           tau = NULL, par.alasso = 1, par.scad = 3.7,
                           thr.delta = 0.01){
  
  result <- switch(method,
                   tapering = cce.taper(cmat, taper, group),
                   glasso = rglasso(cmat, lambda, weight, nlambda, ratio, N),
                   eigen.cleaning = eigen.cleaning(cmat, N),
                   thresholding = cce.thr(cmat, thr.type, thr, tau, par.alasso, par.scad, thr.delta))
  
  return(result)
}

## if matrix inversion is failed, the result is set to NA
myinv <- function(X){
  
  result <- try(solve(X), silent = TRUE)
  
  #if(!is.numeric(result)) result <- try(ginv(X))
  
  if(!is.numeric(result)) result <- NA
  
  return(result)
}

################################## 
#   tapering related functions
##################################

cce.taper <- function(cmat, taper, group){
  
  if(missing(taper)){
    taper <- outer(group, group, FUN = function(x, y) x == y)
  } 
  
  cmat.tilde <- cmat * taper
  
  return(list(cmat = cmat.tilde, pmat = myinv(cmat.tilde)))
}

################################## 
#   glasso related functions
##################################

rglasso.path <- function(S, lambda, weight = TRUE, nlambda, ratio){
  
  if(weight){
    W <- sqrt(diag(S) %o% diag(S))
  }else{
    W <- matrix(1, nrow(S), ncol(S))
  }
  
  diag(W) <- 0
  
  if(missing(lambda)){
    lambda.max <- max(abs((S/W)[upper.tri(S)]))
    lambda.min <- ratio * lambda.max
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  }
  
  f <- function(rho){
    res <- glassoFast::glassoFast(S, rho * W)
    return(list(w = res$w, wi = res$wi))
  }
  
  out <- lapply(lambda, FUN = f)
  w <- lapply(out, function(x) x$w)
  wi <- lapply(out, function(x) x$wi)
  
  return(list(w = w, wi = wi, lambda = lambda))
}

IC.rglasso <- function(out, s, n){
  
  f <- function(wi){
    nll <- -log(det(wi)) + sum(wi * s)
    J <- sum(wi[upper.tri(wi, diag = TRUE)] != 0)
    return(c(aic = nll + (2/n)*J, bic = nll + (log(n)/n)*J))
  }
  
  res <- sapply(out$wi, FUN = f)
  
  return(res)
}

effective.n <- function(yuima, method, frequency, utime){
  
  if(method == "HY" | method == "THY"){
    out <- min(length(yuima))
  }else if(method == "SRC" | method == "SBPC"){
    if(missing(utime)) utime <- set_utime(get.zoo.data(yuima))
    out <- utime/frequency
  }else if(method == "RK" | method == "SIML"){
    out <- min(length(yuima))^(2/5)
  }else if(method == "TSCV"){
    out <- min(length(yuima))^(1/3)
  }else{
    out <- sqrt(min(length(yuima)))
  }
  
  return(out)  
}

rglasso <- function(S, lambda = "aic", weight = TRUE, nlambda = 10, ratio, N){
  
  if(lambda == "aic" | lambda == "bic"){
    
    if(missing(ratio)) ratio <- sqrt(log(nrow(S))/N)
    
    out <- rglasso.path(S, weight = weight, nlambda = nlambda, ratio = ratio)
    ic <- IC.rglasso(out, S, N)
    
    idx <- which.min(ic[lambda, ])
    #w <- out$w[[idx]]
    wi <- out$wi[[idx]]
    w <- solve(wi)
    lambda <- out$lambda[idx]
      
  }else{
    out <- rglasso.path(S, lambda = lambda, weight = weight)
    #w <- out$w[[1]]
    wi <- out$wi[[1]]
    w <- solve(wi)
  }
  
  attr(w, "lambda") <- lambda
  attr(wi, "lambda") <- lambda
  
  dimnames(w) <- dimnames(S)
  dimnames(wi) <- dimnames(S)
  
  return(list(cmat = w, pmat = wi))
}

################################## 
#   eigen.cleaning related functions
##################################

eigen.cleaning <- function(S, N){
  
  p <- ncol(S)
  
  #if(missing(N)) N <- effective.n(yuima, method, frequency, utime)
  
  q <- N/p
  
  R <- cov2cor(S)
  ed <- eigen(R, symmetric = TRUE)
  lambda <- ed$values
  
  sigma2 <- 1 - max(lambda)/p
  lambda.max <- sigma2 * (1 + 1/q + 2/sqrt(q))
  idx <- (lambda > lambda.max)
  delta <- (sum(lambda[lambda > 0]) - sum(lambda[idx]))/(p - sum(idx))
  
  lambda[!idx] <- delta
  
  obj <- sqrt(diag(S))
  cmat <- tcrossprod(ed$vectors %*% diag(lambda), ed$vectors)
  cmat <- obj * cmat %*% diag(obj)
  #pmat <- solve(cmat)
  pmat <- tcrossprod(ed$vectors %*% diag(1/lambda), ed$vectors)
  pmat <- (1/obj) * cmat %*% diag(1/obj)
  
  dimnames(cmat) <- dimnames(S)
  dimnames(pmat) <- dimnames(S)
  
  return(list(cmat = cmat, pmat = pmat))
}

################################## 
#   thresholding related functions
##################################

thr.hard <- function(S, thr){
  out <- S * (abs(S) > thr)
  diag(out) <- diag(S)
  return(out)
}

thr.soft <- function(S, thr){
  out <- sign(S) * pmax(abs(S) - thr, 0)
  diag(out) <- diag(S)
  return(out)
}

thr.alasso <- function(S, thr, par.alasso = 1){
  out <- sign(S) * 
    pmax(abs(S) - thr^(par.alasso + 1)/abs(S)^par.alasso, 0)
  diag(out) <- diag(S)
  return(out)
}  

thr.scad <- function(S, thr, par.scad = 3.7){
  
  out <- S
  
  idx <- (abs(S) <= par.scad * thr)
  out[idx] <- ((par.scad - 1) * S[idx] 
               - sign(S[idx]) * par.scad * thr[idx])/(par.scad - 2)
    
  idx <- (abs(S) <= 2 * thr)
  out[idx] <- sign(S[idx]) * pmax(abs(S[idx]) - thr[idx], 0)
  
  diag(out) <- diag(S)
  
  return(out)
}

thr.fun <- function(S, type = "hard", thr, par.alasso = 1, par.scad = 3.7){
  
  thr <- matrix(thr, nrow(S), ncol(S))
  
  out <- switch (type,
    hard = thr.hard(S, thr),
    soft = thr.soft(S, thr),
    alasso = thr.alasso(S, thr, par.alasso),
    scad = thr.scad(S, thr, par.scad)
  )
  
  return(out)
}

check.pd <- function(S){
  
  lambda <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  
  if(lambda < 0 | isTRUE(all.equal(lambda, 0))){
    return(FALSE)
  }else{
    return(TRUE)
  }
  
}

cce.thr <- function(S, type = "hard", thr = NULL, tau = NULL, par.alasso = 1, par.scad = 3.7, thr.delta = 0.01){
  
  if(is.null(thr)){
    
    if(is.null(tau)){
      
      #if(check.pd(S)){
      #  tau <- 0
      #}else{
      #  
      #  R <- cov2cor(S)
      #  tau.seq <- sort(abs(R[upper.tri(R)])) + .Machine$double.eps
      #  
      #  for(tau in tau.seq){
      #    tmp.R <- thr.fun(R, type, tau, par.alasso, par.scad)
      #    if(check.pd(tmp.R)) break
      #  }
      #}
      R <- cov2cor(S)
      tau.seq <- seq(0, 1, by = thr.delta)
      
      for(tau in tau.seq){
        tmp.R <- thr.fun(R, type, tau, par.alasso, par.scad)
        if(check.pd(tmp.R)) break
      }
      #f <- function(tau){
      #  tmp.R <- thr.fun(R, type, tau, par.alasso, par.scad)
      #  return(min(eigen(tmp.R)$value))
      #}
      #tau <- uniroot(f, interval = c(0, 1))$root
      
    }
    
    thr <- tau * tcrossprod(sqrt(diag(S)))
    
  }
  
  cmat <- thr.fun(S, type, thr, par.alasso, par.scad)
  pmat <- myinv(cmat)
  
  attr(cmat, "thr") <- thr
  attr(pmat, "thr") <- thr
  
  dimnames(cmat) <- dimnames(S)
  dimnames(pmat) <- dimnames(S)
  
  return(list(cmat = cmat, pmat = pmat))
}


