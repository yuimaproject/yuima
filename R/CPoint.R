qmleL <- function(yuima, start, method="BFGS", t, print=FALSE, lower, upper, ...){
	
	call <- match.call()
	
	if( missing(yuima))
	yuima.stop("yuima object is missing.")
	
## param handling
	
## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops	
	if( missing(start) ) 
		yuima.stop("Starting values for the parameters are missing.")
	
	diff.par <- yuima@model@parameter@diffusion
		
	if(!is.list(start))
	yuima.stop("Argument 'start' must be of list type.")
	
	fullcoef <- diff.par
	npar <- length(fullcoef)

	nm <- names(start)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
    start <- start[order(oo)]
    nm <- names(start)
	
	idx.diff <- match(diff.par, nm)
	
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
	
	times <- time(yuima@data@zoo.data[[1]])
	minT <- as.numeric(times[1])
	maxT <- as.numeric(times[length(times)])
	
	if(missing(t))
	 t <- mean(c(minT,maxT))
	k <- max(which(times <=t),na.rm=TRUE)[1]
	
	if(k<2)
	 k <- 2
	if(k>n-2)
	 k <- n-2

	env <- new.env()
#	assign("X",  as.matrix(onezoo(yuima)), env=env)
#	assign("deltaX",  matrix(0, n-1, d.size), env=env)
#	for(t in 1:(n-1))
#	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	assign("X",  as.matrix(onezoo(yuima)[1:k,]), env=env)
	assign("deltaX",  matrix(0, k-1, d.size), env=env)
	for(t in 1:(k-1))
	env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	
	assign("h", deltat(yuima@data@zoo.data[[1]]), env=env)
		
	f <- function(p) {
		mycoef <- as.list(p)
		names(mycoef) <- nm
		sum(pminusquasilogl(yuima=yuima, param=mycoef, print=print, env))
	}
		

	oout <- NULL
	
	
	if(length(start)>1){ # multidimensional optim
		oout <- optim(start, f, method = method, hessian = FALSE, lower=lower, upper=upper)
				
			} else { ### one dimensional optim
				opt <- optimize(f, ...) ## an interval should be provided
                oout <- list(par = opt$minimum, value = opt$objective)
			} ### endif( length(start)>1 )
		

	oout <- oout[c("par","value")]
	return(oout)
}


qmleR <- function(yuima, start, method="BFGS", t, print=FALSE, lower, upper, ...){
	
	call <- match.call()
	
	if( missing(yuima))
	yuima.stop("yuima object is missing.")
	
## param handling
	
## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops	
	if( missing(start) ) 
	yuima.stop("Starting values for the parameters are missing.")
	
	diff.par <- yuima@model@parameter@diffusion
	
	if(!is.list(start))
	yuima.stop("Argument 'start' must be of list type.")
	
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(start)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
    start <- start[order(oo)]
    nm <- names(start)
	
	idx.diff <- match(diff.par, nm)
	
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
	
	times <- time(yuima@data@zoo.data[[1]])
	minT <- as.numeric(times[1])
	maxT <- as.numeric(times[length(times)])
	
	if(missing(t))
	t <- mean(c(minT,maxT))
	k <- max(which(times <=t),na.rm=TRUE)[1]
	
	if(k<2)
	k <- 2
	if(k>n-2)
	k <- n-2


	env <- new.env()

	assign("X",  as.matrix(onezoo(yuima)[-(1:k),]), env=env)
	assign("deltaX",  matrix(0, dim(env$X)[1]-1, d.size), env=env)
	for(t in 1:(dim(env$X)[1]-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	
	assign("h", deltat(yuima@data@zoo.data[[1]]), env=env)
	
	f <- function(p) {
		mycoef <- as.list(p)
		names(mycoef) <- nm
		sum(pminusquasilogl(yuima=yuima, param=mycoef, print=print, env))
	}
	
	
	oout <- NULL
	
	
	if(length(start)>1){ # multidimensional optim
		oout <- optim(start, f, method = method, hessian = FALSE, lower=lower, upper=upper)
		
	} else { ### one dimensional optim
		opt <- optimize(f, ...) ## an interval should be provided
		oout <- list(par = opt$minimum, value = opt$objective)
	} ### endif( length(start)>1 )
	
	
	oout <- oout[c("par","value")]
	return(oout)
}

CPoint <- function(yuima, param1, param2, print=FALSE, plot=FALSE){
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]

	d.size <- yuima@model@equation.number

	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), env=env)
	assign("deltaX",  matrix(0, n-1, d.size), env=env)
	for(t in 1:(n-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), env=env)
	
	QL1 <- pminusquasilogl(yuima=yuima, param=param1, print=print, env)
	QL2 <- pminusquasilogl(yuima=yuima, param=param2, print=print, env)

    D <- sapply(2:(n-1), function(x) sum(QL1[1:x]) + sum(QL2[-(1:x)]))
	D <- c(D[1], D, D[length(D)])
	D <- ts(D, start=0, deltat=deltat(yuima@data@zoo.data[[1]]))
	if(plot)
	 plot(D,type="l", main="change point statistics")
	tau.hat <- index(yuima@data@zoo.data[[1]])[which.min(D)]	

	return(list(tau=tau.hat, param1=param1, param2=param2))
}




# partial quasi-likelihood
# returns a vector of conditionational minus-log-likelihood terms
# the whole negative log likelihood is the sum

pminusquasilogl <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
		yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")

    
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
	n.theta1 <- length(theta1)
	n.theta <- n.theta1
	
	d.size <- yuima@model@equation.number
#n <- length(yuima)[1]
	n <- dim(env$X)[1]
	
	vec <- env$deltaX 
	
	K <- -0.5*d.size * log( (2*pi*h) )
	
	

	QL <- 0
	pn <- numeric(n-1)
	diff <- diffusion.term(yuima, param, env)
	dimB <- dim(diff[, , 1])
	
	if(is.null(dimB)){  # one dimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t]^2
			logdet <- log(yB)
			pn[t] <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
			QL <- QL+pn[t]
			
		}
	} else {  # multidimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t] %*% t(diff[, , t])
			logdet <- log(det(yB))
			if(is.infinite(logdet) ){ # should we return 1e10?
				pn[t] <- log(1)
				yuima.warn("singular diffusion matrix")
				return(1e10)
			}else{
				pn[t] <- K - 0.5*logdet + 
				((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
				QL <- QL+pn[t]
			}
		}
	}
	
	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		QL <- 1e10
	}

	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	
	
	return(-pn)
	
}


pminusquasiloglL <- function(yuima, param, print=FALSE, env){
	
	diff.par <- yuima@model@parameter@diffusion
	fullcoef <- diff.par
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
	
    
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	
	h <- env$h
	
    theta1 <- unlist(param[idx.diff])
	n.theta1 <- length(theta1)
	n.theta <- n.theta1
	
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	vec <- env$deltaX 
	
	K <- -0.5*d.size * log( (2*pi*h) )
	
	
	
	QL <- 0
	pn <- numeric(n-1)
	diff <- diffusion.term(yuima, param, env)
	dimB <- dim(diff[, , 1])
	
	if(is.null(dimB)){  # one dimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t]^2
			logdet <- log(yB)
			pn[t] <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB) 
			QL <- QL+pn[t]
			
		}
	} else {  # multidimensional X
		for(t in 1:(n-1)){
			yB <- diff[, , t] %*% t(diff[, , t])
			logdet <- log(det(yB))
			if(is.infinite(logdet) ){ # should we return 1e10?
				pn[t] <- log(1)
				yuima.warn("singular diffusion matrix")
				return(1e10)
			}else{
				pn[t] <- K - 0.5*logdet + 
				((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]) 
				QL <- QL+pn[t]
			}
		}
	}
	
	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		QL <- 1e10
	}
	
	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	
	return(-pn)
	
}




