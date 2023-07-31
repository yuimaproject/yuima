
## utilities

# function to evaluate a vector of expressions
evalvec <- function(expr, env){
  
  out <- sapply(expr, FUN = eval, envir = env)
  dim(out) <- attr(expr, "dim")
  
  return(out)
}

# function to substitute values for some variables
substvec <- function(s, env){
  
  myfunc <- function(x) 
    do.call("substitute",list(str2lang(x), env)) |> deparse1()
  
  out <- sapply(s, FUN = myfunc)
  dim(out) <- attr(s, "dim")
  
  return(out)
}

# function to check which components are degenerate
is.degenerate <- function(yuima){
  
  diffusion <- yuima@model@diffusion
  d <- length(diffusion)
  out <- logical(d)
  
  for(j in 1:d){
    out[j] <- all(calculus::e2c(diffusion[[j]]) == "(0)")
  }
  
  return(out)
}

# function to extract components of a yuima object
subsetyuima <- function(yuima, idx){
  
  B <- lapply(yuima@model@diffusion[idx], FUN = calculus::e2c)
  B <- do.call("rbind", B)
  
  new.model <- setModel(drift = yuima@model@drift[idx],
                        diffusion = B,
                        hurst = yuima@model@hurst,
                        state.variable = yuima@model@state.variable[idx],
                        time.variable = yuima@model@time.variable,
                        solve.variable = yuima@model@solve.variable[idx],
                        xinit = yuima@model@xinit[idx]
                        )
  
  new.model@noise.number <- yuima@model@noise.number
  
  new.yuima <- setYuima(data = setData(get.zoo.data(yuima)[idx]),
                        model = new.model)
  
  return(new.yuima)
}

# minus quasi-log likelihood function to estimate diffusion parameters
minusH1 <- function(yuima, theta1, idx.x, env){
  
  Cmat <- diffusion.term(yuima, theta1, env)[idx.x, idx.x, ]
  
  h <- env$h
  vec <- env$deltaX
  
  QL <- 0
  
  if(length(idx.x) == 1){  # one dimensional X
    
    for(j in 1:n){
      yB <- Cmat[j]^2
      logdet <- log(yB)
      pn <- -0.5*logdet - 0.5*vec[j]^2/(h*yB)
      QL <- QL + pn
    }
    
  } else {  # multidimensional X
    for(j in 1:n){
      yB <- tcrossprod(Cmat[ , ,j]) 
      logdet <- log(det(yB))
      if(is.infinite(logdet) ){ # should we return 1e10?
        pn <- log(1)
        yuima.warn("singular diffusion matrix")
        return(1e10)
      }else{
        pn <-  - 0.5*logdet +
                 (-1/(2*h)) * crossprod(vec[j, ], solve(yB)) %*% vec[j, ]
        QL <- QL + pn
      }
    }
  }
  
  return(-QL)
}

# minus adaptive quasi-log likelihood function to estimate 
# drift parameters with non-degenerate diffusion terms
minusH2 <- function(yuima, theta2, idx.x, Cmat, env){
  
  h <- env$h
  vec <- env$deltaX
  
  hA <- h * drift.term(yuima, theta2, env)[ ,idx.x]
  
  QL <- 0
  
  if(length(idx.x) == 1){  # one dimensional X
    
    for(j in 1:n){
      yB <- Cmat[j]^2
      pn <- -0.5 * (vec[j] - hA[j])^2/(h*yB)
      QL <- QL + pn
    }
    
  } else {  # multidimensional X
    for(j in 1:n){
      yB <- tcrossprod(Cmat[ , ,j]) 
      pn <-  (-1/(2*h)) * crossprod(vec[j, ] - hA[ ,j], solve(yB)) %*% (vec[j, ] - hA[ ,j])
      QL <- QL + pn
    }
  }
  
  return(-QL)
}


# QMLE for diffusion parameters
init.est.theta1 <- function(yuima, start, method = "L-BFGS-B", #fixed = list(),
                            envir = globalenv(), lower, upper, 
                            idx.x, ...){
  
  call <- match.call() 
  
  diff.par <- yuima@model@parameter@diffusion
  
  nm <- names(start) 
  
  idx.diff <- match(diff.par, nm)
  
  env <- new.env(parent = envir) ##Kurisaki 4/4/2021
  
  Z <- as.matrix(onezoo(yuima))
  dZ <- diff(Z)
  dX <- dZ[ ,idx.x]
  
  assign("X",  as.matrix(onezoo(yuima)), envir = env)
  assign("deltaX",  dX, envir = env)
  assign("h", deltat(get.zoo.data(yuima)[[1]]), envir = env)
  
  f <- function(p) {
    mycoef <- as.list(p)
    names(mycoef) <- nm[idx.diff]
    return(minusH1(yuima, mycoef, idx.x, env))
  }
  
  new.start <- start[idx.diff] # considering only initial guess for diffusion
  names(new.start) <- nm[idx.diff]
  
  mydots <- as.list(call)[-(1:2)]
  mydots$print <- NULL
  mydots$rcpp <- NULL #KK 08/07/16
  mydots$fixed <- NULL
  mydots$fn <- as.name("f")
  mydots$start <- NULL
  mydots$par <- unlist(new.start)
  mydots$hessian <- FALSE
  mydots$upper <- as.numeric(unlist( upper[ nm[idx.diff] ]))
  mydots$lower <- as.numeric(unlist( lower[ nm[idx.diff] ]))
  mydots$joint <- NULL # LM 08/03/16
  mydots$aggregation <- NULL # LM 08/03/16
  mydots$threshold <- NULL #SMI 2/9/14
  mydots$envir <- NULL ##Kurisaki 4/4/2021
  mydots$Est.Incr <- NULL ##Kurisaki 4/10/2021
  mydots$print <- NULL ##Kito 4/17/2021
  mydots$aggregation <- NULL ##Kito 4/17/2021
  mydots$rcpp <- NULL ##Kito 4/17/2021
  mydots$idx.x <- NULL
  
  if((length(mydots$par)>1) | any(is.infinite(c(mydots$upper,mydots$lower)))){
    mydots$method<-method     ##song
    oout <- do.call(optim, args=mydots)
  } else {
    mydots$f <- mydots$fn
    mydots$fn <- NULL
    mydots$par <- NULL
    mydots$hessian <- NULL
    mydots$interval <- as.numeric(c(unlist(lower[diff.par]),unlist(upper[diff.par])))
    mydots$lower <- NULL
    mydots$upper <- NULL
    mydots$method<- NULL
    mydots$envir <- NULL ##Kurisaki 4/4/2021
    mydots$Est.Incr <- NULL ##Kurisaki 4/8/2021
    mydots$print <- NULL ##Kito 4/17/2021
    mydots$aggregation <- NULL ##Kito 4/17/2021
    mydots$rcpp <- NULL ##Kito 4/17/2021
    mydots$idx.x <- NULL
    opt1 <- do.call(optimize, args=mydots)  
    theta1 <- opt1$minimum
    names(theta1) <- diff.par
    oout <- list(par = theta1, value = opt1$objective)
  }
  theta1 <- oout$par
  
  return(theta1)
}

# adaptive QMLE for drift parameters with non-degenerate diffusion terms
init.est.theta2 <- function(yuima, start, method = "L-BFGS-B", #fixed = list(),
                            print = FALSE, envir = globalenv(), 
                            lower, upper, idx.x, theta1.hat, ...){
  
  call <- match.call()
  
  drift.par <- unique(all.vars(yuima@model@drift[idx.x]))
  drift.par <- setdiff(drift.par, c(yuima@model@solve.variable,
                                    yuima@model@state.variable))
  
  #drift.par <- yuima@model@parameter@drift[idx.theta2]
  
  nm <- names(start)
  
  idx.drift <- match(drift.par, nm)
  
  env <- new.env(parent = envir) ##Kurisaki 4/4/2021
  
  ### this part is temporary and should be fixed ####
  drift.par.full <- yuima@model@parameter@drift
  nm.full <- names(drift.par.full)
  for(i in seq_along(drift.par.full)){
    assign(drift.par.full[i], 0, env)
    #assign(drift.par.full[i], start[nm.full[i]], env)
  }
  ###
  
  Z <- as.matrix(onezoo(yuima))
  dZ <- diff(Z)
  dX <- dZ[ ,idx.x]
  
  assign("X",  as.matrix(onezoo(yuima)), envir = env)
  assign("deltaX",  dX, envir = env)
  assign("h", deltat(get.zoo.data(yuima)[[1]]), envir = env)
  
  Cmat <- diffusion.term(yuima, theta1.hat, env)[idx.x, idx.x, ]
  
  f <- function(p) {
    mycoef <- as.list(p)
    names(mycoef) <- nm[idx.drift]
    return(minusH2(yuima, mycoef, idx.x, Cmat, env))
  }
  
  new.start <- start[idx.drift] # considering only initial guess for diffusion
  names(new.start) <- nm[idx.drift]
  
  mydots <- as.list(call)[-(1:2)]
  mydots$print <- NULL
  mydots$rcpp <- NULL #KK 08/07/16
  mydots$fixed <- NULL
  mydots$fn <- as.name("f")
  mydots$start <- NULL
  mydots$par <- unlist(new.start)
  mydots$hessian <- FALSE
  mydots$upper <- as.numeric(unlist( upper[ nm[idx.drift] ]))
  mydots$lower <- as.numeric(unlist( lower[ nm[idx.drift] ]))
  mydots$joint <- NULL # LM 08/03/16
  mydots$aggregation <- NULL # LM 08/03/16
  mydots$threshold <- NULL #SMI 2/9/14
  mydots$envir <- NULL ##Kurisaki 4/4/2021
  mydots$Est.Incr <- NULL ##Kurisaki 4/10/2021
  mydots$print <- NULL ##Kito 4/17/2021
  mydots$aggregation <- NULL ##Kito 4/17/2021
  mydots$rcpp <- NULL ##Kito 4/17/2021
  mydots$idx.x <- NULL
  mydots$theta1.hat <- NULL
  
  if((length(mydots$par)>1) | any(is.infinite(c(mydots$upper,mydots$lower)))){
    mydots$method<-method     ##song
    oout <- do.call(optim, args=mydots)
  } else {
    mydots$f <- mydots$fn
    mydots$fn <- NULL
    mydots$par <- NULL
    mydots$hessian <- NULL
    mydots$interval <- as.numeric(c(unlist(lower[drift.par]),unlist(upper[drift.par])))
    mydots$lower <- NULL
    mydots$upper <- NULL
    mydots$method<- NULL
    mydots$envir <- NULL ##Kurisaki 4/4/2021
    mydots$Est.Incr <- NULL ##Kurisaki 4/8/2021
    mydots$print <- NULL ##Kito 4/17/2021
    mydots$aggregation <- NULL ##Kito 4/17/2021
    mydots$rcpp <- NULL ##Kito 4/17/2021
    mydots$idx.x <- NULL
    mydots$theta1.hat <- NULL
    opt1 <- do.call(optimize, args=mydots)  
    theta2 <- opt1$minimum
    names(theta2) <- drift.par
    oout <- list(par = theta2, value = opt1$objective)
  }
  theta2 <- oout$par
  
  return(theta2)
}


## main part
qmle.degenerate <- function(yuima, start, method = "L-BFGS-B", 
                            #fixed = list(), 
                            envir = globalenv(), 
                            lower, upper, ...){
  
  #### Determine degenerate terms ####
  idx.x <- which(!is.degenerate(yuima))
  if(length(idx.x) == 0) {
    return(qmle(yuima, start = start, method = method, 
                lower = lower, envir = envir, lower = lower, 
                upper = upper, joint = FALSE, ...))
  }
  
  #### Initial estimation ####
  
  theta1.hat <- init.est.theta1(yuima, start = start, 
                                method = method, envir = envir, 
                                lower = lower, upper = upper, 
                                idx.x = idx.x, ...)
  theta2.hat <- init.est.theta2(yuima, start = start,
                                method = method, envir = envir, 
                                lower = lower, upper = upper, 
                                idx.x = idx.x, theta1.hat = theta1.hat,
                                ...)
  
  #### Initial estimation - finish ####
  
  #### Environment setting ####
  
  env <- new.env(parent = envir) ##Kurisaki 4/4/2021
  
  theta1 <- yuima@model@parameter@diffusion
  theta2 <- unique(all.vars(yuima@model@drift[idx.x]))
  theta2 <- setdiff(theta2, c(yuima@model@solve.variable,
                              yuima@model@state.variable))
  
  for(i in seq_along(theta1)){
    assign(theta1[i], theta1.hat[theta1[i]], env)
  }
  
  for(i in seq_along(theta2)){
    assign(theta2[i], theta2.hat[theta2[i]], env)
  }
  
  #assign(theta2, theta2.hat, env)
  
  #assign("h", deltat(get.zoo.data(yuima)[[1]]), envir = env)
  
  #### Environment setting - finish ####
  
  #### Computation of constituents of quasi likelihood  ####
  
  Z <- as.matrix(onezoo(yuima))
  n <- nrow(Z) - 1
  
  dZ <- diff(Z)
  dX <- as.matrix(dZ[ ,idx.x])
  dY <- as.matrix(dZ[ ,-idx.x])
  
  A <- calculus::e2c(yuima@model@drift[idx.x])
  H <- calculus::e2c(yuima@model@drift[-idx.x])
  
  #B <- as.character(unlist(yuima@model@diffusion[idx.x]))
  #B <- matrix(B, nrow = length(idx.x), byrow = TRUE)
  B <- lapply(yuima@model@diffusion[idx.x], FUN = calculus::e2c)
  B <- do.call("rbind", B)
  
  Cmat <- B %mx% t(B)
  
  h <- deltat(get.zoo.data(yuima)[[1]])
  
  svar <- yuima@model@state.variable
  Hx <- calculus::gradient(H, var = svar[idx.x])
  Hy <- calculus::gradient(H, var = svar[-idx.x])
  Hxx <- calculus::hessian(H, var = svar[idx.x])
  
  L.H <- (Hx %mx% A) %sum% 
    (0.5 %prod% (Hxx %dot% Cmat)) %sum%
    (Hy %mx% H)
  
  G <- H %sum% ((0.5 * h) %prod% L.H)
  
  tHx <- t(Hx)
  V <- Hx %mx% Cmat %mx% tHx
  
  Cinv <- calculus::mxinv(Cmat)
  Vinv <- calculus::mxinv(V)
  
  Sinv11 <- Cinv %sum% (3 %prod% tHx %mx% Vinv %mx% Hx)
  Sinv12 <- -6 %prod% tHx %mx% Vinv
  Sinv21 <- t(Sinv12)
  Sinv22 <- 12 %prod% Vinv
  
  Sinv <- rbind(cbind(Sinv11, Sinv12),
                cbind(Sinv21, Sinv22))
  
  detSinv <- calculus::mxdet(Sinv)
  logdetS <- paste0("-log", calculus::wrap(detSinv))
  
  Dj <- c(h^(-1/2) %prod% ("deltaX" %sum% ((-h) %prod% A)),
          h^(-3/2) %prod% ("deltaY" %sum% ((-h) %prod% G)))
  
  H3summand <- 0.5 %prod% ((Dj %mx% Sinv %mx% Dj) %sum% logdetS)
  
  #A <- parse(text = A)
  #G <- parse(text = G)
  #Sinv <- parse(text = Sinv)
  #detSinv <- parse(text = detSinv)
  #A <- calculus::c2e(A)
  #G <- calculus::c2e(G)
  #Sinv <- calculus::c2e(Sinv)
  #detSinv <- calculus::c2e(detSinv)
  eH3summand <- calculus::c2e(H3summand)
  
  #### Computation of constituents of quasi likelihood - finish  ####
  
  #### Identification of unknown parameters to estimate ####
  
  #par <- yuima@model@parameter@drift[idx.theta3]
  par <- unique(all.vars(yuima@model@drift[-idx.x]))
  par <- setdiff(par, c(yuima@model@solve.variable,
                        yuima@model@state.variable))
  
  nm <- names(start)
  
  idx <- match(par, nm)
  new.start <- start[idx] # considering only initial guess for diffusion
  names(new.start) <- nm[idx]
  
  #### Identification of unknown parameters to estimate - finish ####
  
  #### Definition of objective function ####
  f <- function(p){
    
    for(i in 1:length(p)){
      assign(nm[idx[i]], p[i], env)
    }
    #names(p) <- nm[idx]
    
    QL <- 0
    
    for(j in 1:n){
      
      for(i in 1:length(svar)){
        assign(svar[i], Z[j,i], env)
      }
      
      assign("deltaX", dX[j, ], envir = env)
      assign("deltaY", dY[j, ], envir = env)
      
      QL <- QL + eval(eH3summand, env)
      
      #tmpA <- sapply(A, FUN = eval, envir = env)
      #tmpG <- sapply(G, FUN = eval, envir = env)
      #tmpSinv <- sapply(Sinv, FUN = eval, envir = env)
      #tmpSinv <- matrix(tmpSinv, d.size)
      #tmpSinv <- array(tmpSinv, dim = attr(Sinv, "dim"))
      #tmpA <- evaluate(A, var = p, params = as.list(env))
      #tmpG <- evaluate(G, var = p, params = as.list(env))
      #tmpSinv <- evaluate(Sinv, var = p, params = as.list(env))
      #tmpdetSinv <- evaluate(detSinv, var = p, params = as.list(env))
      #tmpA <- evalvec(A, env)
      #tmpG <- evalvec(G, env)
      #tmpSinv <- evalvec(Sinv, env)
      
      #Dj <- c(h^(-1/2) * (dX[j, ] - h * tmpA),
      #        h^(-3/2) * (dY[j, ] - h * tmpG))
      
      #QL <- QL + 0.5 * (Dj %*% tmpSinv %*% Dj - 
      #                    log(eval(detSinv, env)))
      #QL <- QL + 0.5 * (Dj %*% tmpSinv %*% Dj - log(tmpdetSinv))
      
    }
    
    return(QL)
  }
  
  #### Definition of objective function - finish ####
  
  call <- match.call()
  mydots <- as.list(call)[-(1:2)]
  mydots$print <- NULL
  mydots$rcpp <- NULL #KK 08/07/16
  mydots$fixed <- NULL
  mydots$fn <- as.name("f")
  mydots$start <- NULL
  mydots$par <- unlist(new.start)
  mydots$hessian <- FALSE
  mydots$upper <- as.numeric(unlist( upper[ nm[idx] ]))
  mydots$lower <- as.numeric(unlist( lower[ nm[idx] ]))
  mydots$joint <- NULL # LM 08/03/16
  mydots$aggregation <- NULL # LM 08/03/16
  mydots$threshold <- NULL #SMI 2/9/14
  mydots$envir <- NULL ##Kurisaki 4/4/2021
  mydots$Est.Incr <- NULL ##Kurisaki 4/10/2021
  mydots$print <- NULL ##Kito 4/17/2021
  mydots$aggregation <- NULL ##Kito 4/17/2021
  mydots$rcpp <- NULL ##Kito 4/17/2021
  mydots$idx.x <- NULL
  
  if((length(mydots$par)>1) | any(is.infinite(c(mydots$upper,mydots$lower)))){
    mydots$method<-method     ##song
  }else{
    mydots$method <- "Brent"
  }
  
  oout <- do.call(optim, args=mydots)
  theta3 <- oout$par
  names(theta3) <- par
  
  coef <- c(theta1.hat, theta2.hat, theta3)
  
  #### One-Step Estimation ####
  
  ## Compute the gradient and hessian of the summand of H1
  dH1summand <- 
    calculus::gradient(H3summand, var = theta1) |>
    calculus::c2e()
  ddH1summand <- 
    calculus::hessian(H3summand, var = theta1) |>
    calculus::c2e()
  
  ## Compute the summand of H23
  env.onestep <- new.env()
  
  #for(i in seq_along(theta1)){
  #  assign(theta1[i], as.numeric(theta1.hat[theta1[i]]), envir = env.onestep)
  #}
  
  #for(i in seq_along(theta2)){
  #  assign(theta2[i], as.numeric(theta2.hat[theta2[i]]), envir = env.onestep)
  #}
  
  for(i in seq_along(par)){
    assign(par[i], as.numeric(theta3[i]), envir = env.onestep)
  }
  
  Sinv.onestep <- substvec(Sinv, env.onestep)
  #print(Sinv.onestep[1,1])
  H23summand <- 0.5 %prod% Dj %mx% Sinv.onestep %mx% Dj
  
  ## Compute the gradient and hessian of the summand of H23
  dH23summand <- 
    calculus::gradient(H23summand, var = c(theta2, par)) |>
    calculus::c2e()
  ddH23summand <- 
    calculus::hessian(H23summand, var = c(theta2, par)) |>
    calculus::c2e()
  
  ## Evaluation
  dH1 <- double(length(theta1))
  ddH1 <- diag(length(theta1))
  dH23 <- double(length(theta2) + length(par))
  ddH23 <- diag(length(theta2) + length(par))
  
  for(j in 1:n){
    
    for(i in 1:length(svar)){
      assign(svar[i], Z[j,i], env)
    }
    
    assign("deltaX", dX[j, ], envir = env)
    assign("deltaY", dY[j, ], envir = env)
    
    dH1 <- dH1 + evalvec(dH1summand, env)
    ddH1 <- ddH1 + evalvec(ddH1summand, env)
    dH23 <- dH23 + evalvec(dH23summand, env)
    ddH23 <- ddH23 + evalvec(ddH23summand, env)
    
  }
  
  theta1.onestep <- theta1.hat - solve(as.matrix(ddH1), dH1)
  theta23.onestep <- c(theta2.hat, theta3) - 
    solve(as.matrix(ddH23), dH23)
  
  coef.onestep <- c(theta1.onestep, theta23.onestep)
  
  #### One-Step Estimation - finish ####
  
  #### Standard Error Estimation ####
  idx22 <- 1:length(theta2)
  
  ddH2 <- ddH23[idx22,idx22]
  ddH3 <- ddH23[-idx22,-idx22]
  
  #vcov <- bdiag(solve(as.matrix(ddH1)),
  #              solve(as.matrix(ddH23)[idx22,idx22]),
  #              solve(as.matrix(ddH23)[-idx22,-idx22]))
  vcov <- bdiag(solve(as.matrix(ddH1)),
                solve(as.matrix(ddH2)),
                solve(as.matrix(ddH3)))
  #se <- sqrt(diag(vcov))
  
  myfixed <- as.numeric(rep(NA, length(coef.onestep)))
  names(myfixed) <- names(coef.onestep)
  
  final_res<-new("yuima.qmle", 
                 call = call, 
                 coef = coef.onestep, 
                 #fullcoef = unlist(mycoef),
                 fullcoef = coef.onestep,
                 fixed = myfixed,
                 vcov = as.matrix(vcov), 
                 min = oout$value, 
                 details = oout, 
                 minuslogl = f,
                 nobs = as.integer(max(length(yuima)) - 1), 
                 method = method, 
                 model = yuima@model)
  
  #return(list(initial = coef, onestep = coef.onestep,
  #            vcov = vcov, se = se))
  return(final_res)
}

