##::quasi-likelihood function

##::extract drift term from yuima
##::para: parameter of drift term (theta2)

### TO BE FIXED: all caculations should be made on a private environment to
### avoid problems.
### I have rewritten drift.term and diff.term instead of calc.drift and
### calc.diffusion to make them independent of the specification of the 
### parameters.  S.M.I. 22/06/2010


drift.term <- function(yuima, theta, env){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
	DRIFT <- yuima@model@drift
#	n <- length(yuima)[1]
	n <- dim(env$X)[1]
	
	drift <- matrix(0,n,d.size)
	tmp.env <- new.env()
	assign(yuima@model@time.variable, env$time, envir=tmp.env)

	
	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]], envir=tmp.env)
	}
	
	for(d in 1:d.size){
		assign(modelstate[d], env$X[,d], envir=tmp.env)
	}
	for(d in 1:d.size){
		drift[,d] <- eval(DRIFT[d], envir=tmp.env)
	}

	return(drift)  
}


diffusion.term <- function(yuima, theta, env){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
	DIFFUSION <- yuima@model@diffusion
#	n <- length(yuima)[1]
	n <- dim(env$X)[1]
    tmp.env <- new.env()
	assign(yuima@model@time.variable, env$time, envir=tmp.env)
	diff <- array(0, dim=c(d.size, r.size, n))
	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]],envir=tmp.env)
	}

	for(d in 1:d.size){
		assign(modelstate[d], env$X[,d], envir=tmp.env)
	}

	for(r in 1:r.size){
		for(d in 1:d.size){
			diff[d, r, ] <- eval(DIFFUSION[[d]][r], envir=tmp.env)
		}
	}
	return(diff)
}


### I have rewritten qmle as a version of ml.ql
### This function has an interface more similar to mle.
### ml.ql is limited in that it uses fixed names for drift and diffusion
### parameters, while yuima model allows for any names.
### also, I am using the same interface of optim to specify upper and lower bounds
### S.M.I. 22/06/2010


qmle <- function(yuima, start, method="BFGS", fixed = list(), print=FALSE, 
 lower, upper, joint=FALSE, ...){

  if(is(yuima@model, "yuima.carma")&& length(yuima@model@info@scale.par)!=0){
    method<-"L-BFGS-B"
  }
	call <- match.call()
	
	if( missing(yuima))
		yuima.stop("yuima object is missing.")
	
## param handling
	
## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops	
	if( missing(start) ) 
	 yuima.stop("Starting values for the parameters are missing.")

  #14/12/2013 We modify the QMLE function when the model is a Carma(p,q).
  # In this case we use a two step procedure:
  # First) The Coefficient are obtained by QMLE computed using the Kalman Filter.
  # Second) Using the result in Brockwell, Davis and Yang (2007) we retrieve 
  # the underlying Levy. The estimated increments are used to find the Lévy parameters.

#   if(is(yuima@model, "yuima.carma")){
#     yuima.warm("two step procedure for carma(p,q)")
#     return(null)
#   }
#   
  
	diff.par <- yuima@model@parameter@diffusion
	
#	24/12
  if(is(yuima@model, "yuima.carma") && length(diff.par)==0
	   && length(yuima@model@parameter@jump)!=0){
    diff.par<-yuima@model@parameter@jump
	}
  
  if(is(yuima@model, "yuima.carma") && length(yuima@model@parameter@jump)!=0){
    yuima.warn("carma(p,q): the qmle for a carma(p,q) driven by a Jump process will be implemented as soon as possible ")
    return(null)
  }
  
  # 24/12
  if(is(yuima@model, "yuima.carma") && length(yuima@model@info@lin.par)>0){
    yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
    return(null)    
  }
    
  
  drift.par <- yuima@model@parameter@drift
	jump.par <- yuima@model@parameter@jump
	measure.par <- yuima@model@parameter@measure
	common.par <- yuima@model@parameter@common
	
	JointOptim <- joint
	if(length(common.par)>0){
		JointOptim <- TRUE
		yuima.warn("Drift and diffusion parameters must be different. Doing
					  joint estimation, asymptotic theory may not hold true.")
	# 24/12
#     if(is(yuima@model, "yuima.carma")){
#            JointOptim <- TRUE
# 		       yuima.warm("Carma(p.q): The case of common parameters in Drift and Diffusion Term will be implemented as soon as possible,")
# 		       #return(null)
# 		     }
	}

	 if(!is(yuima@model, "yuima.carma")){
    	if(length(jump.par)+length(measure.par)>0)
    		yuima.stop("Cannot estimate the jump models, yet")
	 }
  
  
	if(!is.list(start))
		yuima.stop("Argument 'start' must be of list type.")

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
	if (is(yuima@model, "yuima.carma")){
	  # 24/12
    d.size <-1
	}
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
  if (is(yuima@model, "yuima.carma")){
    #24/12 If we consider a carma model,
    # the observations are only the first column of env$X
	  env$X<-as.matrix(env$X[,1])
	  env$deltaX<-as.matrix(env$deltaX[,1])
	}
  assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 
	
  for(t in 1:(n-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)

	f <- function(p) {
        mycoef <- as.list(p)
#		print(nm[-idx.fixed])
#		print(nm)
		if(length(idx.fixed)>0)
		 names(mycoef) <- nm[-idx.fixed]
		else
		 names(mycoef) <- nm
        mycoef[fixed.par] <- fixed
	    minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }
		
	 fj <- function(p) {
		 mycoef <- as.list(p)
		 names(mycoef) <- nm
		 mycoef[fixed.par] <- fixed
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }

	 oout <- NULL
   HESS <- matrix(0, length(nm), length(nm))
	 colnames(HESS) <- nm
	 rownames(HESS) <- nm
	 HaveDriftHess <- FALSE
	 HaveDiffHess <- FALSE
	 if(length(start)){
		if(JointOptim){ ### joint optimization
			if(length(start)>1){ #Â?multidimensional optim
				oout <- optim(start, fj, method = method, hessian = TRUE, lower=lower, upper=upper)
				HESS <- oout$hessian
				if(is(yuima@model,"yuima.carma") && length(yuima@model@info@scale.par)!=0){
				  b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
				  idx.b0<-match(b0,rownames(HESS))
				  HESS<-HESS[-idx.b0,]
				  HESS<-HESS[,-idx.b0]
				}
				HaveDriftHess <- TRUE
				HaveDiffHess <- TRUE
			} else { ### one dimensional optim
				opt1 <- optimize(f, ...) ## an interval should be provided
                oout <- list(par = opt1$minimum, value = opt1$objective)
			} ### endif( length(start)>1 )
		} else {  ### first diffusion, then drift
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
			mydots$print <- NULL
			mydots$fixed <- NULL
			mydots$fn <- as.name("f")
			mydots$start <- NULL
			mydots$par <- unlist(new.start)
			mydots$hessian <- FALSE
			mydots$upper <- unlist( upper[ nm[idx.diff] ])
			mydots$lower <- unlist( lower[ nm[idx.diff] ])
			if(length(mydots$par)>1){
			 oout <- do.call(optim, args=mydots)
			} else {
			 mydots$f <- mydots$fn
			 mydots$fn <- NULL
			 mydots$par <- NULL
			 mydots$hessian <- NULL	
			 mydots$method <- NULL	
			 mydots$interval <- as.numeric(c(unlist(lower[diff.par]),unlist(upper[diff.par]))) 
			 mydots$lower <- NULL	
			 mydots$upper <- NULL	
			 opt1 <- do.call(optimize, args=mydots)
			 theta1 <- opt1$minimum
			 names(theta1) <- diff.par
			 oout <- list(par = theta1, value = opt1$objective) 
			}
			theta1 <- oout$par
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
			mydots$print <- NULL
			mydots$fixed <- NULL
			mydots$fn <- as.name("f")
			mydots$start <- NULL
			mydots$par <- unlist(new.start)
			mydots$hessian <- FALSE
			mydots$upper <- unlist( upper[ nm[idx.drift] ])
			mydots$lower <- unlist( lower[ nm[idx.drift] ])
			if(length(mydots$par)>1){
			  if(is(yuima@model, "yuima.carma")&& length(yuima@model@info@scale.par)!=0){
			    name_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
          index_b0<-match(name_b0,nm)
			    mydots$lower[index_b0]<-1
          mydots$upper[index_b0]<-1+10^(-7)
			  }
			  oout1 <- do.call(optim, args=mydots)
	#		  oout1 <- optim(mydots$par,f,method = "L-BFGS-B" , lower = mydots$lower, upper = mydots$upper)
			} else {
				mydots$f <- mydots$fn
				mydots$fn <- NULL
				mydots$par <- NULL
				mydots$hessian <- NULL	
				mydots$method <- NULL	
				mydots$interval <- as.numeric(c(lower[drift.par],upper[drift.par])) 
				opt1 <- do.call(optimize, args=mydots)
				theta2 <- opt1$minimum
				names(theta2) <- drift.par
				oout1 <- list(par = theta2, value = as.numeric(opt1$objective)) 	
			}
			theta2 <- oout1$par
			} ## endif(length(idx.drift)>0)
			oout1 <- list(par=  c(theta1, theta2))
			names(oout1$par) <- c(diff.par,drift.par)
			oout <- oout1

		} ### endif JointOptim
    } else {
		list(par = numeric(0L), value = f(start))
	}
	
	 	
	 fDrift <- function(p) {
		 mycoef <- as.list(p)
		 names(mycoef) <- drift.par
		 mycoef[diff.par] <- coef[diff.par]
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }

	 fDiff <- function(p) {
		 mycoef <- as.list(p)
		 names(mycoef) <- diff.par
		 mycoef[drift.par] <- coef[drift.par]
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }
	 
	 coef <- oout$par
	 control=list()
	 par <- coef
  if(!is(yuima@model,"yuima.carma")){
	  names(par) <- c(diff.par, drift.par)
       nm <- c(diff.par, drift.par)
  }
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
     if(is(yuima@model,"yuima.carma") && length(yuima@model@info@scale.par)!=0){
       b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
       idx.b0<-match(b0,rownames(HESS))
       HESS<-HESS[-idx.b0,]
       HESS<-HESS[,-idx.b0]
     }
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
	
	min <- minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
  
  dummycov<-matrix(0,length(coef),length(coef))
  rownames(dummycov)<-names(coef)
  colnames(dummycov)<-names(coef)
  dummycov[rownames(vcov),colnames(vcov)]<-vcov
  vcov<-dummycov
  
    new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
       vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl, 
       method = method)
}


quasilogl <- function(yuima, param, print=FALSE){

	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(yuima:::onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	for(t in 1:(n-1))
	env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 

	-minusquasilogl(yuima=yuima, param=param, print=print, env)
}


minusquasilogl <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	#  24/12
	if(is(yuima@model, "yuima.carma") && length(yuima@model@info@lin.par)==0
	   && length(yuima@model@parameter@jump)!=0){
	  diff.par<-yuima@model@parameter@jump
	}
	
	# 24/12
	if(is(yuima@model, "yuima.carma") && length(yuima@model@info@lin.par)>0  ){
	  yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
	  return(null)    
	}
	
  
  drift.par <- yuima@model@parameter@drift
	
	fullcoef <- NULL
	
	if(length(diff.par)>0)
	fullcoef <- diff.par
	
	if(length(drift.par)>0)
	fullcoef <- c(fullcoef, drift.par)
	
#	fullcoef <- c(diff.par, drift.par)
	npar <- length(fullcoef)
#	cat("\nparam\n")
#	print(param)
#	cat("\nfullcoef\n")
#	print(fullcoef)
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
    theta2 <- unlist(param[idx.drift])
	n.theta1 <- length(theta1)
	n.theta2 <- length(theta2)
	n.theta <- n.theta1+n.theta2
	
	d.size <- yuima@model@equation.number
	
  
	n <- length(yuima)[1]
	
  
  if (is(yuima@model, "yuima.carma")){
	  # 24/12
	  d.size <-1
    # We build the two step procedure as described in
    if(length(yuima@model@info@scale.par)!=0){
       prova<-as.numeric(param)
       names(prova)<-fullcoef[oo]
       param<-prova[c(length(prova):1)]
       
       y<-as.numeric(env$X)
       tt<-env$time
       p<-yuima@model@info@p
       a<-param[c(1:p)]
#        a_names<-names(param[c(1:p)])
#        names(a)<-a_names
       b<-param[c((p+1):(length(param)-p+1))]
#        b_names<-names(param[c((p+1):(length(param)-p+1))])
#        names(b)<-b_names
       if(length(b)==1){
         b<-1
       } else{ 
         indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
         b[indx_b0]<-1
       }
       sigma<-tail(param,1)
       strLog<-yuima.carma.loglik1(y, tt, a, b, sigma)
       QL<-strLog$loglikCdiag
      }else {
        yuima.warn("carma(p,q): the scale parameter is equal to 1. We will implemented as soon as possible")
        return(NULL)
    }
	} else{   
  	drift <- drift.term(yuima, param, env)
  	diff <- diffusion.term(yuima, param, env)
  	
  	QL <- 0
  	
  	pn <- 0
  
  	
  	vec <- env$deltaX-h*drift[-n,]
  
  	K <- -0.5*d.size * log( (2*pi*h) )
  
  	dimB <- dim(diff[, , 1])
  
  	if(is.null(dimB)){  # one dimensional X
  	  for(t in 1:(n-1)){
  		yB <- diff[, , t]^2
  		logdet <- log(yB)
  		pn <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
  		QL <- QL+pn
  			
  		}
  	} else {  # multidimensional X
  	 for(t in 1:(n-1)){
  		yB <- diff[, , t] %*% t(diff[, , t])
  		logdet <- log(det(yB))
  		if(is.infinite(logdet) ){ # should we return 1e10?
  			pn <- log(1)
  			yuima.warn("singular diffusion matrix")
  			return(1e10)
  		}else{
  			pn <- K - 0.5*logdet + 
  					  ((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
  			QL <- QL+pn
  		}
  	 }
  	}
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




MatrixA<-function (a) 
{
  #Build Matrix A in the state space representation of Carma(p,q) 
  #given the autoregressive coefficient
  pp = length(a)
  af = cbind(rep(0, pp - 1), diag(pp - 1))
  af = rbind(af, -a[pp:1])
  return(af)
}


yuima.Vinfinity<-function(elForVInf,v){
  # We find the infinity stationary variance-covariance matrix
  A<-elForVInf$A
  sigma<-elForVInf$sigma
  p<-dim(A)[1]  
  matrixV<-matrix(0,p,p)
  matrixV[upper.tri(matrixV,diag=TRUE)]<-v
  matrixV[lower.tri(matrixV)]<-matrixV[upper.tri(matrixV)]
  l<-rbind(matrix(rep(0,p-1),p-1,1),1)
  
  RigSid<-l%*%t(l)
  Matrixobj<-A%*%matrixV+matrixV%*%t(A)+sigma^2*RigSid
  obj<-sum(Matrixobj^2)
  obj
}


carma.kalman<-function(y, tt, p, q, a,bvector, sigma){
  
  V_inf0<-matrix(diag(rep(1,p)),p,p)
  
  A<-MatrixA(a)
  u<-diff(tt)[1]
  expA<-as.matrix(expm(A*u))
  v<-as.numeric(V_inf0[upper.tri(V_inf0,diag=TRUE)])
  elForVInf<-list(A=A,sigma=sigma)
  
  V_inf_vect<-nlm(yuima.Vinfinity, v, elForVInf = elForVInf)$estimate
  V_inf<-matrix(0,p,p)
  V_inf[upper.tri(V_inf,diag=TRUE)]<-V_inf_vect
  V_inf[lower.tri(V_inf)]<-V_inf[upper.tri(V_inf)]
  
  V_inf[abs(V_inf)<=1.e-06]=0
  
  
  SIGMA_err<-V_inf-expA%*%V_inf%*%t(expA)
  
  statevar0<-matrix(rep(0, p),p,1)
  Qmatr<-SIGMA_err
  
  # set
  statevar<-statevar0
  
  #   SigMatr<-expA%*%V_inf%*%t(expA)+Qmatr
  
  SigMatr<-Qmatr
  #SigMatr<-V_inf
  
  zc<-matrix(bvector,1,p)
  loglstar = 0
  loglstar1 = 0
  for(t in 1:length(tt)){ 
    # prediction
    statevar<-expA%*%statevar
    SigMatr<-expA%*%SigMatr%*%t(expA)+Qmatr
    # forecast
    Uobs<-y[t]-zc%*%statevar
    sd_2<-zc%*%SigMatr%*%t(zc)
    
    #correction
    Kgain<-SigMatr%*%t(zc)%*%solve(sd_2)
    #  SigMatr<-SigMatr-Kgain%*%as.numeric(sd_2)%*%t(Kgain)
    
    statevar<-statevar+Kgain%*%Uobs
    SigMatr<-SigMatr-Kgain%*%zc%*%SigMatr
    # construction of log-likelihood
    #     loglstar1<-loglstar1+log(dnorm(as.numeric(Uobs), mean = 0, sd = sqrt(as.numeric(sd_2))))
    #     sdsig<-sqrt(as.numeric(sd_2))
    term_int<--0.5*(log((2*pi)^(length(Uobs))*det(sd_2))+t(Uobs)%*%solve(sd_2)%*%Uobs)
    loglstar<-loglstar+term_int
  }
  return(list(loglstar=as.numeric(loglstar),s2hat=as.numeric(sd_2)))
}



yuima.carma.loglik1<-function (y, tt, a, b, sigma) 
{
  #This code compute the LogLik using kalman filter
  
  # if(a_0!=0){we need to correct the Y_t for the mean}
  # if(sigma!=1){we need to write}
  p <- as.integer(length(a))
  
  bvector <- rep(0, p)
  q <- length(b)
  bvector <- c(b, rep(0, p - q))
  
  
  sigma<-sigma
  y<-y
  
  xxalt<-carma.kalman(y, tt, p, q, a,bvector,sigma)
  list(loglikCdiag = xxalt$loglstar,s2hat=xxalt$s2hat)
}


loglik5 <- function(param) {
  a<-param[1:pp]
  b <- 1
  if(qq>0){
    b<-param[(pp + 1):(pp + qq)]
  }
  
  sigma<-tail(param,n=1)
  xx <- yuima.carma.loglik1(y, tt, a, b, sigma)
  
  return(xx$loglikCdiag)
}








# returns the vector of log-transitions instead of the final quasilog
quasiloglvec <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	
	fullcoef <- NULL
	
	if(length(diff.par)>0)
	fullcoef <- diff.par
	
	if(length(drift.par)>0)
	fullcoef <- c(fullcoef, drift.par)
	
	npar <- length(fullcoef)

	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
    theta2 <- unlist(param[idx.drift])
	n.theta1 <- length(theta1)
	n.theta2 <- length(theta2)
	n.theta <- n.theta1+n.theta2
	
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	
	drift <- drift.term(yuima, param, env)
	diff <- diffusion.term(yuima, param, env)
	
	QL <- numeric(n-1)  ## here is the difference
	
	pn <- 0
    
	
	vec <- env$deltaX-h*drift[-n,]
    
	K <- -0.5*d.size * log( (2*pi*h) )
    
	dimB <- dim(diff[, , 1])
    
	if(is.null(dimB)){  # one dimensional X
        for(t in 1:(n-1)){
            yB <- diff[, , t]^2
            logdet <- log(yB)
            pn <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
            QL[t] <- pn
			
		}
	} else {  # multidimensional X
        for(t in 1:(n-1)){
            yB <- diff[, , t] %*% t(diff[, , t])
            logdet <- log(det(yB))
            if(is.infinite(logdet) ){ # should we return 1e10?
                pn <- log(1)
                yuima.warn("singular diffusion matrix")
                return(1e10)
            }else{
                pn <- K - 0.5*logdet + 
                ((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
                QL[t] <- pn
            }
        }
	}
	return(QL)
}




## test function using nlmnb instead of optim. Not mush difference

qmle2 <- function(yuima, start, method="BFGS", fixed = list(), print=FALSE, 
lower, upper, joint=FALSE, ...){
    
	call <- match.call()
	
	if( missing(yuima))
    yuima.stop("yuima object is missing.")
	
    ## param handling
	
    ## FIXME: maybe we should choose initial values at random within lower/upper
    ##        at present, qmle stops	
	if( missing(start) ) 
    yuima.stop("Starting values for the parameters are missing.")
    
	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	jump.par <- yuima@model@parameter@jump
	measure.par <- yuima@model@parameter@measure
	common.par <- yuima@model@parameter@common
	
	JointOptim <- joint
	if(length(common.par)>0){
		JointOptim <- TRUE
		yuima.warn("Drift and diffusion parameters must be different. Doing
        joint estimation, asymptotic theory may not hold true.")
	}
    
    
	if(length(jump.par)+length(measure.par)>0)
    yuima.stop("Cannot estimate the jump models, yet")
	
	if(!is.list(start))
    yuima.stop("Argument 'start' must be of list type.")
    
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
	
	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 
	for(t in 1:(n-1))
    env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
    
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
    
	f <- function(p) {
        mycoef <- as.list(p)
        #		print(nm[-idx.fixed])
        #		print(nm)
		if(length(idx.fixed)>0)
        names(mycoef) <- nm[-idx.fixed]
		else
        names(mycoef) <- nm
        mycoef[fixed.par] <- fixed
	    minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }
    
    fj <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- nm
        mycoef[fixed.par] <- fixed
        minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }
    
    oout <- NULL
    HESS <- matrix(0, length(nm), length(nm))
    colnames(HESS) <- nm
    rownames(HESS) <- nm
    HaveDriftHess <- FALSE
    HaveDiffHess <- FALSE
    if(length(start)){
		if(JointOptim){ ### joint optimization
			if(length(start)>1){ #Â?multidimensional optim
				oout <- nlminb(start, fj, lower=lower, upper=upper)
                print(oout)
				HESS <- oout$hessian
				HaveDriftHess <- TRUE
				HaveDiffHess <- TRUE
			} else { ### one dimensional optim
				opt1 <- optimize(f, ...) ## an interval should be provided
                oout <- list(par = opt1$minimum, value = opt1$objective)
			} ### endif( length(start)>1 )
		} else {  ### first diffusion, then drift
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
                mydots$print <- NULL
                mydots$fixed <- NULL
                mydots$fn <- as.name("f")
                mydots$objective <- mydots$fn
                mydots$start <- NULL
                mydots$par <- unlist(new.start)
                mydots$start <- unlist(new.start)
                mydots$hessian <- NULL
                mydots$upper <- unlist( upper[ nm[idx.diff] ])
                mydots$lower <- unlist( lower[ nm[idx.diff] ])
                if(length(mydots$par)>1){
                    mydots$fn <- NULL
                    mydots$par <- NULL

                    oout <- do.call(nlminb, args=mydots)
                    #                    oout <- do.call(optim, args=mydots)
                } else {
                    mydots$f <- mydots$fn
                    mydots$fn <- NULL
                    mydots$par <- NULL
                    mydots$hessian <- NULL	
                    mydots$method <- NULL	
                    mydots$start <- NULL 
                    mydots$objective <- NULL
                    mydots$interval <- as.numeric(c(unlist(lower[diff.par]),unlist(upper[diff.par]))) 
                    mydots$lower <- NULL	
                    mydots$upper <- NULL	
                    opt1 <- do.call(optimize, args=mydots)
                    theta1 <- opt1$minimum
                    names(theta1) <- diff.par
                    oout <- list(par = theta1, value = opt1$objective) 
                }
                theta1 <- oout$par
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
                mydots$print <- NULL
                mydots$fixed <- NULL
                mydots$fn <- as.name("f")
                mydots$objective <- mydots$fn
                mydots$start <- NULL
                mydots$par <- unlist(new.start)
                mydots$start <- unlist(new.start)
                mydots$hessian <- NULL
                mydots$upper <- unlist( upper[ nm[idx.drift] ])
                mydots$lower <- unlist( lower[ nm[idx.drift] ])
                
                if(length(mydots$par)>1){
                    mydots$par <- NULL
                    

mydots$fn <- NULL
                    oout1 <- do.call(nlminb, args=mydots)
                    #                    oout1 <- do.call(optim, args=mydots)
                } else {
                    mydots$f <- mydots$fn
                    mydots$fn <- NULL
                    mydots$par <- NULL
                    mydots$start <- NULL
                    mydots$objective <- NULL
                    mydots$hessian <- NULL	
                    mydots$method <- NULL	
                    mydots$interval <- as.numeric(c(lower[drift.par],upper[drift.par])) 
                    opt1 <- do.call(optimize, args=mydots)
                    theta2 <- opt1$minimum
                    names(theta2) <- drift.par
                    oout1 <- list(par = theta2, value = as.numeric(opt1$objective)) 	
                }
                theta2 <- oout1$par
			} ## endif(length(idx.drift)>0)
			oout1 <- list(par=  c(theta1, theta2))
			names(oout1$par) <- c(diff.par,drift.par)
			oout <- oout1
            
		} ### endif JointOptim
    } else {
		list(par = numeric(0L), value = f(start))
	}
	
    
    fDrift <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- drift.par
        mycoef[diff.par] <- coef[diff.par]
        minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }
    
    fDiff <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- diff.par
        mycoef[drift.par] <- coef[drift.par]
        minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
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
	
	min <- minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	
    new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
    vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl, 
    method = method)
}


