##::quasi-bayes function

#new minusquasilogl "W1","W2" like lse function.

setGeneric("lseBayes",
function(yuima, start,prior,lower,upper, method="mcmc",mcmc=1000,rate=1.0,algorithm="randomwalk")
standardGeneric("lseBayes")
)
setMethod("lseBayes", "yuima",
function(yuima, start,prior,lower,upper, method="mcmc",mcmc=1000,rate=1.0,algorithm="randomwalk")
{
 
  rcpp <- TRUE
  
  joint <- FALSE
  fixed <- numeric(0)
  print <- FALSE
  
  call <- match.call()
  
  if( missing(yuima))
    yuima.stop("yuima object is missing.")
  
  ## param handling
  
  ## FIXME: maybe we should choose initial values at random within lower/upper
  ##        at present, qmle stops	
  
  if(missing(lower) || missing(upper)){
    yuima.stop("lower or upper is missing.")
  }
  
  diff.par <- yuima@model@parameter@diffusion
  drift.par <- yuima@model@parameter@drift
  jump.par <- yuima@model@parameter@jump
  measure.par <- yuima@model@parameter@measure
  common.par <- yuima@model@parameter@common
  
  ## BEGIN Prior construction
  if(!missing(prior)){
    priorLower = numeric(0)
    priorUpper = numeric(0)
    pdlist <- numeric(length(yuima@model@parameter@all))
    names(pdlist) <- yuima@model@parameter@all
    for(i in 1: length(pdlist)){
      if(prior[[names(pdlist)[i]]]$measure.type=="code"){
        expr <- prior[[names(pdlist)[i]]]$df
        code <- suppressWarnings(sub("^(.+?)\\(.+", "\\1", expr, perl=TRUE))
        args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", expr, perl=TRUE)), ","))
        pdlist[i] <- switch(code,
                            dunif=paste("function(z){return(dunif(z, ", args[2], ", ", args[3],"))}"),
                            dnorm=paste("function(z){return(dnorm(z,", args[2], ", ", args[3], "))}"),
                            dbeta=paste("function(z){return(dbeta(z, ", args[2], ", ", args[3], "))}"),
                            dgamma=paste("function(z){return(dgamma(z, ", args[2], ", ", args[3], "))}"),
                            dexp=paste("function(z){return(dexp(z, ", args[2], "))}")
        )
        qf <- switch(code,
                     dunif=paste("function(z){return(qunif(z, ", args[2], ", ", args[3],"))}"),
                     dnorm=paste("function(z){return(qnorm(z,", args[2], ", ", args[3], "))}"),
                     dbeta=paste("function(z){return(qbeta(z, ", args[2], ", ", args[3], "))}"),
                     dgamma=paste("function(z){return(qgamma(z, ", args[2], ", ", args[3], "))}"),
                     dexp=paste("function(z){return(qexp(z, ", args[2], "))}")
        )
        priorLower = append(priorLower,eval(parse("text"=qf))(0.00))
        priorUpper = append(priorUpper,eval(parse("text"=qf))(1.00))
        
        
      }
      
    }
    if(sum(unlist(lower)<priorLower) + sum(unlist(upper)>priorUpper) > 0){
      yuima.stop("lower&upper of prior are out of parameter space.")
    }
    
    names(lower) <- names(pdlist)
    names(upper) <- names(pdlist)
    
    
    
    pd <- function(param){
      value <- 1
      for(i in 1:length(pdlist)){
        value <- value*eval(parse(text=pdlist[[i]]))(param[[i]])
      }
      return(value)
    }
  }else{
    pd <- function(param) return(1)
  }
  ## END Prior construction
  
  
  JointOptim <- joint
  if(length(common.par)>0){
    JointOptim <- TRUE
    yuima.warn("Drift and diffusion parameters must be different. Doing
               joint estimation, asymptotic theory may not hold true.")
  }
  
  
  if(length(jump.par)+length(measure.par)>0)
    yuima.stop("Cannot estimate the jump models, yet")
  
  
  fullcoef <- NULL
  
  if(length(diff.par)>0)
    fullcoef <- diff.par
  
  if(length(drift.par)>0)
    fullcoef <- c(fullcoef, drift.par)
  
  npar <- length(fullcoef)
  
  fixed.par <- names(fixed)
  
  if (any(!(fixed.par %in% fullcoef))) 
    yuima.stop("Some named arguments in 'fixed' are not arguments to the supplied yuima model")
  
  nm <- names(start)
  oo <- match(nm, fullcoef)
  if(any(is.na(oo))) 
    yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
  start <- start[order(oo)]
  if(!missing(prior)){
    pdlist <- pdlist[order(oo)]
  }
  nm <- names(start)
  
  idx.diff <- match(diff.par, nm)
  idx.drift <- match(drift.par, nm)
  idx.fixed <- match(fixed.par, nm)
  tmplower <- as.list( rep( -Inf, length(nm)))
  names(tmplower) <- nm	
  if(!missing(lower)){
    idx <- match(names(lower), names(tmplower))
    if(any(is.na(idx)))
      yuima.stop("names in 'lower' do not match names fo parameters")
    tmplower[ idx ] <- lower	
  }
  lower <- tmplower
  
  tmpupper <- as.list( rep( Inf, length(nm)))
  names(tmpupper) <- nm	
  if(!missing(upper)){
    idx <- match(names(upper), names(tmpupper))
    if(any(is.na(idx)))
      yuima.stop("names in 'lower' do not match names fo parameters")
    tmpupper[ idx ] <- upper	
  }
  upper <- tmpupper
  
  
  
  
  d.size <- yuima@model@equation.number
  n <- length(yuima)[1]
  
  G <- rate
  if(G<=0 || G>1){
    yuima.stop("rate G should be 0 < G <= 1")
  }
  n_0 <- floor(n^G)
  if(n_0 < 2) n_0 <- 2
  #######data is reduced to n_0 before qmle(16/11/2016)
  env <- new.env()
  #assign("X",  yuima@data@original.data[1:n_0,], envir=env)
  assign("X",  as.matrix(onezoo(yuima)[1:n_0,]), envir=env)
  assign("deltaX",  matrix(0, n_0 - 1, d.size), envir=env)
  assign("crossdx",matrix(0,n_0 - 1,d.size*d.size),envir=env) ####(deltaX)%*%t(deltaX).this is used in W1.
  assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env)
  
  assign("Cn.r", rep(1,n_0 - 1), envir=env)
  
  for(t in 1:(n_0 - 1)){
    env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
    env$crossdx[t,] <- as.vector(tcrossprod(env$deltaX[t,]))
  }
  
  assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)

  pp<-0 
  while(1){
    if(n*env$h^pp < 0.1) break
    pp <- pp + 1
  }
  qq <- max(pp,2/G) 
  
  C.temper.diff <- n_0^(2/(qq*G)-1) #this is used in pg.
  C.temper.drift <- (n_0*env$h)^(2/(qq*G)-1) #this is used in pg.

  
  mle <- qmle(yuima, "start"=start, "lower"=lower,"upper"=upper, "method"="L-BFGS-B",rcpp=rcpp)
  start <- as.list(mle@coef)
  
 
  integ <- function(idx.fixed=NULL,f=f,start=start,par=NULL,hessian=FALSE,upper,lower){
    if(length(idx.fixed)==0){
      intf <- adaptIntegrate(f,lowerLimit=lower,upperLimit=upper,fDim=(length(upper)+1))$integral
    }else{
      intf <- adaptIntegrate(f,lowerLimit=lower[-idx.fixed],upperLimit=upper[-idx.fixed],fDim=(length(upper[-idx.fixed])+1))$integral
    }
    return(intf[-1]/intf[1])
  }
  mcinteg <- function(idx.fixed=NULL,f=f,p,start=start,par=NULL,hessian=FALSE,upper,lower,mean,vcov,mcmc){
    if(length(idx.fixed)==0){
      intf <- mcIntegrate(f,p,lowerLimit=lower,upperLimit=upper,mean,vcov,mcmc)
    }else{
      intf <- mcIntegrate(f,p,lowerLimit=lower[-idx.fixed],upperLimit=upper[-idx.fixed],mean[-idx.fixed],vcov[-idx.fixed,-idx.fixed],mcmc)
    }
    return(intf)
  }
  
  mcIntegrate <- function(f,p, lowerLimit, upperLimit,mean,vcov,mcmc){
    
    
    if(algorithm=="randomwalk"){
      x_c <- mean
      p_c <- p(mean)
      val <- f(x_c)

      if(length(mean)>1){
        x <- rmvnorm(mcmc-1,mean,vcov)
        q <- dmvnorm(x,mean,vcov)
        q_c <- dmvnorm(mean,mean,vcov) 
      }else{
        x <- rnorm(mcmc-1,mean,sqrt(vcov))
        q <- dnorm(x,mean,sqrt(vcov))
        q_c <- dnorm(mean,mean,sqrt(vcov)) 
      }
      
      for(i in 1:(mcmc-1)){
        if(length(mean)>1){x_n <- x[i,]}else{x_n <- x[i]}
        if(sum(x_n<lowerLimit)==0 & sum(x_n>upperLimit)==0){
          q_n <- q[i]
          p_n <- p(x_n)
          #u <- runif(1)
          #a <- (p_n*q_c)/(p_c*q_n)
          u <- log(runif(1))
          a <- p_n-p_c+log(q_c/q_n)
          if(u<a){
            p_c <- p_n
            q_c <- q_n
            x_c <- x_n
          }
        }
        val <- val+f(x_c)
      }
      return(unlist(val/mcmc))
    }
    else if(tolower(algorithm)=="mpcn"){ #MpCN
      x_n <- mean
      val <- mean
      logLik_old <- p(x_n)+0.5*length(mean)*log(sqnorm(x_n-mean))
      
      for(i in 1:(mcmc-1)){
        prop <- makeprop(mean,x_n,unlist(lowerLimit),unlist(upperLimit))
        logLik_new <- p(prop)+0.5*length(mean)*log(sqnorm(prop-mean))
        u <- log(runif(1))
        if( logLik_new-logLik_old > u){
          nx_ <- prop
          logLik_old <- logLik_new
        }
        val <- val+f(x_n)
      }
      return(unlist(val/mcmc))
    }
  }
  
  
  #print(mle@coef)
  flagNotPosDif <- 0
  for(i in 1:npar){
    if(mle@vcov[i,i] <= 0) flagNotPosDif <- 1 #Check mle@vcov is positive difinite matrix
  }
  if(flagNotPosDif == 1){
    mle@vcov <- diag(c(rep(1 / n_0,length(diff.par)),rep(1 / (n_0 * env$h),length(drift.par)))) # Redifine mle@vcov
  }
  
  tmpW1 <- minusquasilogl_W1(yuima=yuima, param=mle@coef, print=print, env,rcpp=rcpp)
  tmpW2 <- minusquasilogl_W2(yuima=yuima, param=mle@coef, print=print, env,rcpp=rcpp)
  
  g <- function(p,fixed,idx.fixed){
    mycoef <- mle@coef
    if(length(idx.fixed)>0){
      mycoef[-idx.fixed] <- p
      mycoef[idx.fixed] <- fixed
    }else{
      names(mycoef) <- nm
    }

    if(sum(idx.diff==idx.fixed)>0){
      return(c(1,p)*exp(-minusquasilogl_W1(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmpW1)*pd(param=mycoef))
    }else{
      return(c(1,p)*exp(-minusquasilogl_W2(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmpW2)*pd(param=mycoef))
    }
    
  }
  
  pg <- function(p,fixed,idx.fixed){
    mycoef <- start
    if(length(idx.fixed)>0){
      mycoef[-idx.fixed] <- p
      mycoef[idx.fixed] <- fixed
    }else{
      names(mycoef) <- nm
    }
    
    if(sum(idx.diff==idx.fixed)>0){
      return(C.temper.diff*(-minusquasilogl_W1(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmpW1+log(pd(param=mycoef))))#log
    }else{
      return(C.temper.drift*(-minusquasilogl_W2(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmpW2+log(pd(param=mycoef))))#log
    }
  }
  
  idf <- function(p){return(p)}
  
  #	 fj <- function(p) {
  #		 mycoef <- as.list(p)
  #		 names(mycoef) <- nm
  #		 mycoef[fixed.par] <- fixed
  #		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
  #	 }
  
  oout <- NULL
  HESS <- matrix(0, length(nm), length(nm))
  colnames(HESS) <- nm
  rownames(HESS) <- nm
  HaveDriftHess <- FALSE
  HaveDiffHess <- FALSE
  if(length(start)){
    #		if(JointOptim){ ### joint optimization
    #			if(length(start)>1){ #multidimensional optim
    #				oout <- optim(start, fj, method = method, hessian = TRUE, lower=lower, upper=upper)
    #				HESS <- oout$hessian
    #				HaveDriftHess <- TRUE
    #				HaveDiffHess <- TRUE
    #			} else { ### one dimensional optim
    #				opt1 <- optimize(f, ...) ## an interval should be provided
    #				opt1 <- list(par=integ(f=f,upper=upper,lower=lower,fDim=length(lower)+1),objective=0)
    #               oout <- list(par = opt1$minimum, value = opt1$objective)
    #			} ### endif( length(start)>1 )
    #		} else {  ### first diffusion, then drift
    theta1 <- NULL
    
    old.fixed <- fixed 
    old.start <- start
    
    if(length(idx.diff)>0){
      ## DIFFUSION ESTIMATIOn first
      old.fixed <- fixed
      old.start <- start
      new.start <- start[idx.diff] # considering only initial guess for diffusion
      new.fixed <- fixed
      if(length(idx.drift)>0)	
        new.fixed[nm[idx.drift]] <- start[idx.drift]
      fixed <- new.fixed
      fixed.par <- names(fixed)
      idx.fixed <- match(fixed.par, nm)
      names(new.start) <- nm[idx.diff]
      mydots <- as.list(call)[-(1:2)]
      mydots$fixed <- NULL
      mydots$fn <- as.name("f")
      mydots$start <- NULL
      mydots$par <- unlist(new.start)
      mydots$hessian <- FALSE
      mydots$upper <- unlist( upper[ nm[idx.diff] ])
      mydots$lower <- unlist( lower[ nm[idx.diff] ])
      f <- function(p){return(g(p,fixed,idx.fixed))}
      pf <- function(p){return(pg(p,fixed,idx.fixed))}
      if(length(mydots$par)>1){
        #			 oout <- do.call(optim, args=mydots)
        if(method=="mcmc"){
          oout <- list(par=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=diag(diag(mle@vcov)),mcmc=mcmc))
        }else{
          oout <- list(par=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower,start=start))
        }
      } else {
        mydots$f <- mydots$fn
        mydots$fn <- NULL
        mydots$par <- NULL
        mydots$hessian <- NULL	
        mydots$method <- NULL	
        mydots$interval <- as.numeric(c(unlist(lower[diff.par]),unlist(upper[diff.par]))) 
        mydots$lower <- NULL	
        mydots$upper <- NULL	
        #			 opt1 <- do.call(optimize, args=mydots)
        if(method=="mcmc"){
          opt1 <- list(minimum=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=diag(diag(mle@vcov)),mcmc=mcmc))
        }else{
          opt1 <- list(minimum=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
        }
        theta1 <- opt1$minimum
        names(theta1) <- diff.par
        #			 oout <- list(par = theta1, value = opt1$objective) 
        oout <- list(par=theta1,value=0)
      }
      theta1 <- oout$par
      #names(theta1) <- nm[idx.diff]
      names(theta1) <- diff.par
    } ## endif(length(idx.diff)>0)
    
    theta2 <- NULL
    
    if(length(idx.drift)>0){
      ## DRIFT estimation with first state diffusion estimates
      fixed <- old.fixed
      start <- old.start
      new.start <- start[idx.drift] # considering only initial guess for drift
      new.fixed <- fixed
      new.fixed[names(theta1)] <- theta1
      fixed <- new.fixed
      fixed.par <- names(fixed)
      idx.fixed <- match(fixed.par, nm)
      names(new.start) <- nm[idx.drift]
      
      mydots <- as.list(call)[-(1:2)]
      mydots$fixed <- NULL
      mydots$fn <- as.name("f")
      mydots$start <- NULL
      mydots$par <- unlist(new.start)
      mydots$hessian <- FALSE
      mydots$upper <- unlist( upper[ nm[idx.drift] ])
      mydots$lower <- unlist( lower[ nm[idx.drift] ])
      f <- function(p){return(g(p,fixed,idx.fixed))}
      pf <- function(p){return(pg(p,fixed,idx.fixed))}
      
      if(length(mydots$par)>1){
        #			  oout1 <- do.call(optim, args=mydots)
        if(method=="mcmc"){
          oout1 <- list(par=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=diag(diag(mle@vcov)),mcmc=mcmc))
        }else{
          oout1 <- list(par=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
        }
      } else {
        mydots$f <- mydots$fn
        mydots$fn <- NULL
        mydots$par <- NULL
        mydots$hessian <- NULL	
        mydots$method <- NULL	
        mydots$interval <- as.numeric(c(lower[drift.par],upper[drift.par])) 
        #				opt1 <- do.call(optimize, args=mydots)
        if(method=="mcmc"){
          opt1 <- list(minimum=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=diag(diag(mle@vcov)),mcmc=mcmc))
        }else{
          opt1 <- list(minimum=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
        }
        theta2 <- opt1$minimum
        names(theta2) <- drift.par
        oout1 <- list(par = theta2, value = as.numeric(opt1$objective)) 	
      }
      theta2 <- oout1$par
    } ## endif(length(idx.drift)>0)
    oout1 <- list(par=  c(theta1, theta2))
    names(oout1$par) <- c(diff.par,drift.par)
    oout <- oout1
    
    #		} ### endif JointOptim
  } else {
    list(par = numeric(0L), value = f(start))
  }
  
  
  fDrift <- function(p) {
    mycoef <- as.list(p)
    names(mycoef) <- drift.par
    mycoef[diff.par] <- coef[diff.par]
    minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
  }
  
  fDiff <- function(p) {
    mycoef <- as.list(p)
    names(mycoef) <- diff.par
    mycoef[drift.par] <- coef[drift.par]
    minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
  }
  
  coef <- oout$par
  control=list()
  par <- coef
  names(par) <- c(diff.par, drift.par)
  nm <- c(diff.par, drift.par)
  
  #	 print(par)
  #	 print(coef)
  conDrift <- list(trace = 5, fnscale = 1, 
                   parscale = rep.int(5, length(drift.par)), 
                   ndeps = rep.int(0.001, length(drift.par)), maxit = 100L, 
                   abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
                   beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
                   factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  conDiff <- list(trace = 5, fnscale = 1, 
                  parscale = rep.int(5, length(diff.par)), 
                  ndeps = rep.int(0.001, length(diff.par)), maxit = 100L, 
                  abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
                  beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
                  factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  
  #	 nmsC <- names(con)
  #	 if (method == "Nelder-Mead") 
  #	 con$maxit <- 500
  #	 if (method == "SANN") {
  #		 con$maxit <- 10000
  #		 con$REPORT <- 100
  #	 }
  #	 con[(namc <- names(control))] <- control
  #	 if (length(noNms <- namc[!namc %in% nmsC])) 
  #	 warning("unknown names in control: ", paste(noNms, collapse = ", "))
  #	 if (con$trace < 0) 
  #	 warning("read the documentation for 'trace' more carefully")
  #	 else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 
  #			  0) 
  #	 stop("'trace != 0' needs 'REPORT >= 1'")
  #	 if (method == "L-BFGS-B" && any(!is.na(match(c("reltol", 
  #													"abstol"), namc)))) 
  #	 warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
  #	 npar <- length(par)
  #	 if (npar == 1 && method == "Nelder-Mead") 
  #	 warning("one-diml optimization by Nelder-Mead is unreliable: use optimize")
  #	 
  if(!HaveDriftHess & (length(drift.par)>0)){
    #hess2 <- .Internal(optimhess(coef[drift.par], fDrift, NULL, conDrift))
    hess2 <- optimHess(coef[drift.par], fDrift, NULL, control=conDrift)
    HESS[drift.par,drift.par] <- hess2	 
  }
  
  if(!HaveDiffHess  & (length(diff.par)>0)){
    #hess1 <- .Internal(optimhess(coef[diff.par], fDiff, NULL, conDiff))
    hess1 <- optimHess(coef[diff.par], fDiff, NULL, control=conDiff)
    HESS[diff.par,diff.par] <- hess1	 
  }
  
  oout$hessian <- HESS
  
  vcov <- if (length(coef)) 
    solve(oout$hessian)
  else matrix(numeric(0L), 0L, 0L)
  
  mycoef <- as.list(coef)
  names(mycoef) <- nm
  mycoef[fixed.par] <- fixed
  
  min <- minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
  
  new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
      #       vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl, 
      vcov = vcov,  details = oout, 
      method = method)
  }
)
minusquasilogl_W1 <- function(yuima, param, print=FALSE, env,rcpp=T){ #new logl estimates volatility
  
  diff.par <- yuima@model@parameter@diffusion
  
  drift.par <- yuima@model@parameter@drift
  if(0){
    if(length(yuima@model@info@scale.par)!=0){
      xinit.par <- yuima@model@parameter@xinit
    }
  }
  
  
  if(0 && length(yuima@model@info@lin.par)==0
     && length(yuima@model@parameter@jump)!=0){
    diff.par<-yuima@model@parameter@jump
    # measure.par<-yuima@model@parameter@measure
  }
  
  if(0 && length(yuima@model@info@lin.par)==0
     && length(yuima@model@parameter@measure)!=0){
    measure.par<-yuima@model@parameter@measure
  }
  
  # 24/12
  if(0 && length(yuima@model@info@lin.par)>0  ){
    yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
    return(NULL)
  }
  
  if(0){
    xinit.par <- yuima@model@parameter@xinit
  }
  
  
  drift.par <- yuima@model@parameter@drift
  
  fullcoef <- NULL
  
  if(length(diff.par)>0)
    fullcoef <- diff.par
  
  if(length(drift.par)>0)
    fullcoef <- c(fullcoef, drift.par)
  
  if(0){
    if(length(xinit.par)>0)
      fullcoef <- c(fullcoef, xinit.par)
  }
  
  if(0 && (length(yuima@model@parameter@measure)!=0))
    fullcoef<-c(fullcoef, measure.par)
  
  if(0){
    if("mean.noise" %in% names(param)){
      mean.noise<-"mean.noise"
      fullcoef <- c(fullcoef, mean.noise)
      NoNeg.Noise<-TRUE
    }
  }
  
  
  npar <- length(fullcoef)
  
  nm <- names(param)
  oo <- match(nm, fullcoef)
  
  if(any(is.na(oo)))
    yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
  param <- param[order(oo)]
  nm <- names(param)
  
  idx.diff <- match(diff.par, nm)
  idx.drift <- match(drift.par, nm)
  
  
  if(0){
    idx.xinit <-as.integer(na.omit(match(xinit.par, nm)))
  }
  
  h <- env$h
  
  Cn.r <- env$Cn.r
  
  theta1 <- unlist(param[idx.diff])
  theta2 <- unlist(param[idx.drift])
  
  
  n.theta1 <- length(theta1)
  n.theta2 <- length(theta2)
  n.theta <- n.theta1+n.theta2
  
  
  if(0){
    theta3 <- unlist(param[idx.xinit])
    n.theta3 <- length(theta3)
    n.theta <- n.theta1+n.theta2+n.theta3
  }
  
  
  d.size <- yuima@model@equation.number
  
  
  n <- length(yuima)[1]
  
  
  if (0){
    # 24/12
    d.size <-1
    # We build the two step procedure as described in
    #  if(length(yuima@model@info@scale.par)!=0){
    prova<-as.numeric(param)
    #names(prova)<-fullcoef[oo]
    names(prova)<-names(param)
    param<-prova[c(length(prova):1)]
    time.obs<-env$time.obs
    y<-as.numeric(env$X)
    u<-env$h
    p<-env$p
    q<-env$q
    #         p<-yuima@model@info@p
    ar.par <- yuima@model@info@ar.par
    name.ar<-paste0(ar.par, c(1:p))
    # 	  q <- yuima@model@info@q
    ma.par <- yuima@model@info@ma.par
    name.ma<-paste0(ma.par, c(0:q))
    if (length(yuima@model@info@loc.par)==0){
      
      a<-param[name.ar]
      #        a_names<-names(param[c(1:p)])
      #        names(a)<-a_names
      b<-param[name.ma]
      #        b_names<-names(param[c((p+1):(length(param)-p+1))])
      #        names(b)<-b_names
      if(length(yuima@model@info@scale.par)!=0){
        if(length(b)==1){
          b<-1
        } else{
          indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
          b[indx_b0]<-1
        }
        sigma<-tail(param,1)
      }else {sigma<-1}
      NoNeg.Noise<-FALSE
      if(0){
        if("mean.noise" %in% names(param)){
          
          NoNeg.Noise<-TRUE
        }
      }
      if(NoNeg.Noise==TRUE){
        if (length(b)==p){
          #mean.noise<-param[mean.noise]
          # Be useful for carma driven by a no negative levy process
          mean.y<-mean(y)
          #mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
          #param[mean.noise]<-mean.y/(tail(b,n=1)/tail(a,n=1)*sigma)
        }else{
          mean.y<-0
        }
        y<-y-mean.y
      }
      # V_inf0<-matrix(diag(rep(1,p)),p,p)
      V_inf0<-env$V_inf0
      p<-env$p
      q<-env$q
      strLog<-yuima.carma.loglik1(y, u, a, b, sigma,time.obs,V_inf0,p,q)
    }else if (!rcpp){
      # 01/01
      #          ar.par <- yuima@model@info@ar.par
      #          name.ar<-paste0(ar.par, c(1:p))
      a<-param[name.ar]
      #          ma.par <- yuima@model@info@ma.par
      #          q <- yuima@model@info@q
      name.ma<-paste0(ma.par, c(0:q))
      b<-param[name.ma]
      if(length(yuima@model@info@scale.par)!=0){
        if(length(b)==1){
          b<-1
        } else{
          indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
          b[indx_b0]<-1
        }
        scale.par <- yuima@model@info@scale.par
        sigma <- param[scale.par]
      } else{sigma <- 1}
      loc.par <- yuima@model@info@loc.par
      mu <- param[loc.par]
      
      NoNeg.Noise<-FALSE
      if(0){
        if("mean.noise" %in% names(param)){
          
          NoNeg.Noise<-TRUE
        }
      }
      
      # Lines 883:840 work if we have a no negative noise
      if(0&&(NoNeg.Noise==TRUE)){
        if (length(b)==p){
          mean.noise<-param[mean.noise]
          # Be useful for carma driven by levy process
          #   mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
          mean.y<-mean(y-mu)
          
        }else{
          mean.y<-0
        }
        y<-y-mean.y
      }
      
      
      y.start <- y-mu
      #V_inf0<-matrix(diag(rep(1,p)),p,p)
      V_inf0<-env$V_inf0
      p<-env$p
      q<-env$q
      strLog<-yuima.carma.loglik1(y.start, u, a, b, sigma,time.obs,V_inf0,p,q)
    }
    
    QL<-strLog$loglikCdiag
    #       }else {
    #         yuima.warn("carma(p,q): the scale parameter is equal to 1. We will implemented as soon as possible")
    #         return(NULL)
    #     }
  } else {
    drift_name <- yuima@model@drift
    diffusion_name <- yuima@model@diffusion
    ####data <- yuima@data@original.data
    data <- env$X
    
    thetadim <- length(yuima@model@parameter@all)
    
    noise_number <- yuima@model@noise.number
    
    assign(yuima@model@time.variable,env$time[-length(env$time)])
    for(i in 1:d.size) assign(yuima@model@state.variable[i], data[-length(data[,1]),i])
    for(i in 1:thetadim) assign(names(param)[i], param[[i]])
    
    d_b <- NULL
    for(i in 1:d.size){
      if(length(eval(drift_name[[i]]))==(length(data[,1])-1)){
        d_b[[i]] <- drift_name[[i]] #this part of model includes "x"(state.variable)
      }
      else{
        if(is.na(c(drift_name[[i]][2]))){ #ex. yuima@model@drift=expression(0) (we hope "expression((0))")
          drift_name[[i]] <- parse(text=paste(sprintf("(%s)", drift_name[[i]])))[[1]]
        }
        d_b[[i]] <- parse(text=paste("(",drift_name[[i]][2],")*rep(1,length(data[,1])-1)",sep=""))
        #vectorization
      }
    }
    
    v_a<-matrix(list(NULL),d.size,noise_number)
    for(i in 1:d.size){
      for(j in 1:noise_number){
        if(length(eval(diffusion_name[[i]][[j]]))==(length(data[,1])-1)){
          v_a[[i,j]] <- diffusion_name[[i]][[j]] #this part of model includes "x"(state.variable)
        }
        else{
          if(is.na(c(diffusion_name[[i]][[j]][2]))){
            diffusion_name[[i]][[j]] <- parse(text=paste(sprintf("(%s)", diffusion_name[[i]][[j]])))[[1]]
          }
          v_a[[i,j]] <- parse(text=paste("(",diffusion_name[[i]][[j]][2],")*rep(1,length(data[,1])-1)",sep=""))
          #vectorization
        }
      }
    }
    
    dx_set <- as.matrix((data-rbind(numeric(d.size),as.matrix(data[-length(data[,1]),])))[-1,])
    crossdx_set <- env$crossdx
    
    drift_set <- diffusion_set <- NULL
    #for(i in 1:thetadim) assign(names(param)[i], param[[i]])
    for(i in 1:d.size) drift_set <- cbind(drift_set,eval(d_b[[i]]))
    for(i in 1:noise_number){
      for(j in 1:d.size) diffusion_set <- cbind(diffusion_set,eval(v_a[[j,i]]))
    }
    QL <- W1(crossdx_set,drift_set,diffusion_set,env$h)*(-0.5*env$h*env$h)
  }
  
  
  if(!is.finite(QL)){
    yuima.warn("quasi likelihood is too small to calculate.")
    return(1e10)
  }
  if(print==TRUE){
    yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
  }
  if(is.infinite(QL)) return(1e10)
  return(as.numeric(-QL))
  
}

minusquasilogl_W2 <- function(yuima, param, print=FALSE, env,rcpp=T){#new logl estimates drift
  
  diff.par <- yuima@model@parameter@diffusion
  
  drift.par <- yuima@model@parameter@drift
  if(0){
    if(length(yuima@model@info@scale.par)!=0){
      xinit.par <- yuima@model@parameter@xinit
    }
  }
  
  
  if(0 && length(yuima@model@info@lin.par)==0
     && length(yuima@model@parameter@jump)!=0){
    diff.par<-yuima@model@parameter@jump
    # measure.par<-yuima@model@parameter@measure
  }
  
  if(0 && length(yuima@model@info@lin.par)==0
     && length(yuima@model@parameter@measure)!=0){
    measure.par<-yuima@model@parameter@measure
  }
  
  # 24/12
  if(0 && length(yuima@model@info@lin.par)>0  ){
    yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
    return(NULL)
  }
  
  if(0){
    xinit.par <- yuima@model@parameter@xinit
  }
  
  
  drift.par <- yuima@model@parameter@drift
  
  fullcoef <- NULL
  
  if(length(diff.par)>0)
    fullcoef <- diff.par
  
  if(length(drift.par)>0)
    fullcoef <- c(fullcoef, drift.par)
  
  if(0){
    if(length(xinit.par)>0)
      fullcoef <- c(fullcoef, xinit.par)
  }
  
  if(0 && (length(yuima@model@parameter@measure)!=0))
    fullcoef<-c(fullcoef, measure.par)
  
  if(0){
    if("mean.noise" %in% names(param)){
      mean.noise<-"mean.noise"
      fullcoef <- c(fullcoef, mean.noise)
      NoNeg.Noise<-TRUE
    }
  }
  
  
  npar <- length(fullcoef)
  
  nm <- names(param)
  oo <- match(nm, fullcoef)
  
  if(any(is.na(oo)))
    yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
  param <- param[order(oo)]
  nm <- names(param)
  
  idx.diff <- match(diff.par, nm)
  idx.drift <- match(drift.par, nm)
  
  
  if(0){
    idx.xinit <-as.integer(na.omit(match(xinit.par, nm)))
  }
  
  h <- env$h
  
  Cn.r <- env$Cn.r
  
  theta1 <- unlist(param[idx.diff])
  theta2 <- unlist(param[idx.drift])
  
  
  n.theta1 <- length(theta1)
  n.theta2 <- length(theta2)
  n.theta <- n.theta1+n.theta2
  
  
  if(0){
    theta3 <- unlist(param[idx.xinit])
    n.theta3 <- length(theta3)
    n.theta <- n.theta1+n.theta2+n.theta3
  }
  
  
  d.size <- yuima@model@equation.number
  
  
  n <- length(yuima)[1]
  
  
  if (0){
    # 24/12
    d.size <-1
    # We build the two step procedure as described in
    #  if(length(yuima@model@info@scale.par)!=0){
    prova<-as.numeric(param)
    #names(prova)<-fullcoef[oo]
    names(prova)<-names(param)
    param<-prova[c(length(prova):1)]
    time.obs<-env$time.obs
    y<-as.numeric(env$X)
    u<-env$h
    p<-env$p
    q<-env$q
    #         p<-yuima@model@info@p
    ar.par <- yuima@model@info@ar.par
    name.ar<-paste0(ar.par, c(1:p))
    # 	  q <- yuima@model@info@q
    ma.par <- yuima@model@info@ma.par
    name.ma<-paste0(ma.par, c(0:q))
    if (length(yuima@model@info@loc.par)==0){
      
      a<-param[name.ar]
      #        a_names<-names(param[c(1:p)])
      #        names(a)<-a_names
      b<-param[name.ma]
      #        b_names<-names(param[c((p+1):(length(param)-p+1))])
      #        names(b)<-b_names
      if(length(yuima@model@info@scale.par)!=0){
        if(length(b)==1){
          b<-1
        } else{
          indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
          b[indx_b0]<-1
        }
        sigma<-tail(param,1)
      }else {sigma<-1}
      NoNeg.Noise<-FALSE
      if(0){
        if("mean.noise" %in% names(param)){
          
          NoNeg.Noise<-TRUE
        }
      }
      if(NoNeg.Noise==TRUE){
        if (length(b)==p){
          #mean.noise<-param[mean.noise]
          # Be useful for carma driven by a no negative levy process
          mean.y<-mean(y)
          #mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
          #param[mean.noise]<-mean.y/(tail(b,n=1)/tail(a,n=1)*sigma)
        }else{
          mean.y<-0
        }
        y<-y-mean.y
      }
      # V_inf0<-matrix(diag(rep(1,p)),p,p)
      V_inf0<-env$V_inf0
      p<-env$p
      q<-env$q
      strLog<-yuima.carma.loglik1(y, u, a, b, sigma,time.obs,V_inf0,p,q)
    }else if (!rcpp){
      # 01/01
      #          ar.par <- yuima@model@info@ar.par
      #          name.ar<-paste0(ar.par, c(1:p))
      a<-param[name.ar]
      #          ma.par <- yuima@model@info@ma.par
      #          q <- yuima@model@info@q
      name.ma<-paste0(ma.par, c(0:q))
      b<-param[name.ma]
      if(length(yuima@model@info@scale.par)!=0){
        if(length(b)==1){
          b<-1
        } else{
          indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
          b[indx_b0]<-1
        }
        scale.par <- yuima@model@info@scale.par
        sigma <- param[scale.par]
      } else{sigma <- 1}
      loc.par <- yuima@model@info@loc.par
      mu <- param[loc.par]
      
      NoNeg.Noise<-FALSE
      if(0){
        if("mean.noise" %in% names(param)){
          
          NoNeg.Noise<-TRUE
        }
      }
      
      # Lines 883:840 work if we have a no negative noise
      if(0&&(NoNeg.Noise==TRUE)){
        if (length(b)==p){
          mean.noise<-param[mean.noise]
          # Be useful for carma driven by levy process
          #   mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
          mean.y<-mean(y-mu)
          
        }else{
          mean.y<-0
        }
        y<-y-mean.y
      }
      
      
      y.start <- y-mu
      #V_inf0<-matrix(diag(rep(1,p)),p,p)
      V_inf0<-env$V_inf0
      p<-env$p
      q<-env$q
      strLog<-yuima.carma.loglik1(y.start, u, a, b, sigma,time.obs,V_inf0,p,q)
    }
    
    QL<-strLog$loglikCdiag
    #       }else {
    #         yuima.warn("carma(p,q): the scale parameter is equal to 1. We will implemented as soon as possible")
    #         return(NULL)
    #     }
  } else {
    drift_name <- yuima@model@drift
    diffusion_name <- yuima@model@diffusion
    ####data <- yuima@data@original.data
    data <- env$X
    
    thetadim <- length(yuima@model@parameter@all)
    
    noise_number <- yuima@model@noise.number
    
    assign(yuima@model@time.variable,env$time[-length(env$time)])
    for(i in 1:d.size) assign(yuima@model@state.variable[i], data[-length(data[,1]),i])
    for(i in 1:thetadim) assign(names(param)[i], param[[i]])
    
    d_b <- NULL
    for(i in 1:d.size){
      if(length(eval(drift_name[[i]]))==(length(data[,1])-1)){
        d_b[[i]] <- drift_name[[i]] #this part of model includes "x"(state.variable)
      }
      else{
        if(is.na(c(drift_name[[i]][2]))){ #ex. yuima@model@drift=expression(0) (we hope "expression((0))")
          drift_name[[i]] <- parse(text=paste(sprintf("(%s)", drift_name[[i]])))[[1]]
        }
        d_b[[i]] <- parse(text=paste("(",drift_name[[i]][2],")*rep(1,length(data[,1])-1)",sep=""))
        #vectorization
      }
    }
    
    dx_set <- as.matrix((data-rbind(numeric(d.size),as.matrix(data[-length(data[,1]),])))[-1,])
    drift_set <- diffusion_set <- NULL
    
    for(i in 1:d.size) drift_set <- cbind(drift_set,eval(d_b[[i]]))
    
    QL <- W2(dx_set,drift_set,env$h)*(-0.5*env$h)
  }
  
  
  if(!is.finite(QL)){
    yuima.warn("quasi likelihood is too small to calculate.")
    return(1e10)
  }
  if(print==TRUE){
    yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
  }
  if(is.infinite(QL)) return(1e10)
  return(as.numeric(-QL))
  
}
