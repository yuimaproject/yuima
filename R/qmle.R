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
#	modelpara <- yuima@model@parameter@drift
	DRIFT <- yuima@model@drift
	n <- length(yuima)[1]
	drift <- matrix(0,n,d.size)
#	X <- as.matrix(onezoo(yuima))

	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]])
	}
	for(t in 1:n){
#		Xt <- X[t,]
		for(d in 1:d.size){
#			assign(modelstate[d],Xt[d])
			assign(modelstate[d], env$X[t,d])
		}
		for(d in 1:d.size){
			drift[t,d] <- eval(DRIFT[d])
		}
	}
	return(drift)  
}



diffusion.term <- function(yuima, theta, env){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
	DIFFUSION <- yuima@model@diffusion
	n <- length(yuima)[1]
	diff <- array(0, dim=c(d.size, r.size, n))
#	X <- as.matrix(onezoo(yuima))
	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]])
	}

	for(r in 1:r.size){
		for(t in 1:n){
#			Xt <- X[t, ]
			for(d in 1:d.size){
				assign(modelstate[d], env$X[t,d])
			}
			for(d in 1:d.size){
				diff[d, r, t] <- eval(DIFFUSION[[d]][r])
			}
		}
	}
	return(diff)
}




##::calculate diffusion%*%t(diffusion) matrix
calc.B <- function(diff){
  d.size <- dim(diff)[1]
  n <- dim(diff)[3]
  B <- array(0, dim=c(d.size, d.size, n))
  for(t in 1:n){
    B[, , t] <- diff[, , t]%*%t(diff[, , t])
  }
  return(B)
}




### I have rewritten qmle as a version of ml.ql
### This function has an interface more similar to mle.
### ml.ql is limited in that it uses fixed names for drift and diffusion
### parameters, while yuima model allows for any names.
### also, I am using the same interface of optim to specify upper and lower bounds
### S.M.I. 22/06/2010


qmle <- function(yuima, start, method="BFGS", fixed = list(), print=FALSE, ...){

	call <- match.call()
	
	if( missing(yuima))
		yuima.stop("yuima object is missing.")
	
## param handling
	if( missing(start) )
	 yuima.warn("Starting values for the parameters are missing. Using random initial values.")

	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	jump.par <- yuima@model@parameter@jump
	measure.par <- yuima@model@parameter@measure
	common.par <- yuima@model@parameter@common
	
	if(length(common.par)>0)
		yuima.stop("Drift and diffusion parameters must be different.")
	
	if(length(jump.par)+length(measure.par)>0)
		yuima.stop("Cannot estimate the jump models, yet")
	
	if(!is.list(start))
		yuima.stop("Argument 'start' must be of list type.")

	fullcoef <- c(diff.par, drift.par)
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

	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(yuima:::onezoo(yuima)), env=env)
	assign("deltaX",  matrix(0, n-1, d.size), env=env)
	for(t in 1:(n-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	assign("h", deltat(yuima@data@zoo.data[[1]]), env=env)

	
	f <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- nm
        mycoef[fixed.par] <- fixed
        minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }
		
	oout <- if(length(start)){ 
		optim(start, f, method = method, hessian = TRUE, ...)
    } else {
		list(par = numeric(0L), value = f(start))
	}
	
	coef <- oout$par
    vcov <- if (length(coef)) 
	solve(oout$hessian)
    else matrix(numeric(0L), 0L, 0L)
    min <- oout$value
	
  	mycoef <- as.list(coef)
	names(mycoef) <- nm
	mycoef[fixed.par] <- fixed
	
    new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
       vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl, 
       method = method)
}


quasilogl <- function(yuima, param, print=FALSE){

	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(yuima:::onezoo(yuima)), env=env)
	assign("deltaX",  matrix(0, n-1, d.size), env=env)
	for(t in 1:(n-1))
	env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), env=env)
	
	-minusquasilogl(yuima=yuima, param=param, print=print, env)
}


minusquasilogl <- function(yuima, param, print=FALSE, env){

	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	fullcoef <- c(diff.par, drift.par)
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

	B <- calc.B(diff)
	
	QL <- 0

	pn <- 0
	for(t in 1:(n-1)){
		yB <- as.matrix(B[, , t])
		ydet <- det(yB)
		if(abs(ydet) <1e-7){ # should we return 1e10?
			pn <- log(1)
			yuima.warn("singular diffusion matrix")
			return(1e10)
		}else{
			pn <- log( 1/((2*pi*h)^(d.size/2)*ydet^(1/2)) *
						 exp((-1/(2*h))*t(env$deltaX[t, ]-h*drift[t, ])%*%solve(yB)%*%(env$deltaX[t,]-h*drift[t, ])) )
			QL <- QL+pn
		}
	}
	if(QL==-Inf){
		yuima.warn("quasi likelihood is too small to calculate.")
	}
	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	if(is.infinite(QL)) return(1e10)
	return(-QL)

}



