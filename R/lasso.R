# Initial version of lasso estimation for SDEs
# fixes a small bug in passing arguments to the
# inner optim function

lasso <- function (yuima, lambda0, start, delta = 1, ...) 
{
    call <- match.call()
    if (missing(yuima)) 
    yuima.stop("yuima object 'yuima' is missing.")
    if (missing(lambda0)) {
        pars <- yuima@model@parameter@all
        lambda0 <- rep(1, length(pars))
        names(lambda0) <- pars
        lambda0 <- as.list(lambda0)
    }
    if (missing(start)) 
    yuima.stop("Starting values for the parameters are missing.")
    fail <- lapply(lambda0, function(x) as.numeric(NA))
    cat("\nLooking for MLE estimates...\n")
    fit <- try(qmle(yuima, start = start, ...), silent = TRUE)
    if (class(fit) == "try-error") 
    return(list(mle = fail, sd.mle = NA, lasso = fail, sd.lasso = NA))
    theta.mle <- coef(fit)
    SIGMA <- try(sqrt(diag(vcov(fit))), silent = TRUE)
    if (class(SIGMA) == "try-error") 
    return(list(mle = theta.mle, sd.mle = NA, lasso = fail, sd.lasso = NA))
    H <- try(solve(vcov(fit)), silent = TRUE)
    if (class(H) == "try-error") 
    return(list(mle = theta.mle, sd.mle = SIGMA, lasso = fail, sd.lasso = NA))
    
    lambda <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)^delta
    
    #lambda1 <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)
    idx <- which(lambda > 10000)
    lambda[idx] <- 10000
    f2 <- function(theta) as.numeric(t(theta - theta.mle) %*% 
    H %*% (theta - theta.mle) + lambda %*% abs(theta))
    cat("\nPerforming LASSO estimation...\n")
    args <- list(...)
    args$joint <- NULL
    args$par <- theta.mle
    args$fn <- f2
    args$hessian <- TRUE
    args$delta <- NULL
    args$control <- list(maxit = 30000, temp = 2000, REPORT = 500)
    #  fit2 <- try(optim(theta.mle, f2, hessian = TRUE, args, control = list(maxit = 30000, 
    #     temp = 2000, REPORT = 500)), silent = FALSE)
    fit2 <- try( do.call(optim, args = args), silent = TRUE)
    
    
    if (class(fit2) == "try-error") 
    return(list(mle = theta.mle, sd.mle = SIGMA, lasso = fail, sd.lasso = NA))
    theta.lasso <- fit2$par
    SIGMA1 <- try(sqrt(diag(solve(fit2$hessian))), silent = TRUE)
    if (class(SIGMA1) == "try-error") 
    return(list(mle = theta.mle, sd.mle = SIGMA, lasso = theta.lasso, 
    sd.lasso = NA))
    return(list(mle = theta.mle, sd.mle = SIGMA, lasso = theta.lasso, 
    sd.lasso = SIGMA1, call = call, lambda0 = lambda0))
}

# removed version
old.lasso <- function(yuima, lambda0, start, delta=1, ...){
	
	call <- match.call()
	
	if( missing(yuima))
		yuima.stop("yuima object 'yuima' is missing.")
	
	if( missing(lambda0) ){
	 pars <- yuima@model@parameter@all	
	 lambda0 <- rep(1, length(pars))
	 names(lambda0) <- pars
	 lambda0 <- as.list(lambda0)	
	}
	
## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops	
	if( missing(start) ) 
		yuima.stop("Starting values for the parameters are missing.")
	
	fail <- lapply(lambda0, function(x) as.numeric(NA))
	
	cat("\nLooking for MLE estimates...\n")
	fit <- try(qmle(yuima, start=start,...), silent=TRUE)
	if(class(fit)=="try-error")
	return(list(mle=fail, sd.mle=NA, lasso=fail, sd.lasso=NA))
	
	SIGMA <- try( sqrt(diag(vcov(fit))), silent=TRUE)
	if(class(SIGMA)=="try-error")
	return(list(mle=fail, sd.mle=NA, lasso=fail, sd.lasso=NA))
	
	
	theta.mle <- coef(fit)
	
	H <- try( solve(vcov(fit)), silent=TRUE) 
	
	if(class(H)=="try-error")
	return(list(mle=fail, sd.mle=NA, lasso=fail, sd.lasso=NA))
	
	
#	lambda <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)
	lambda <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)^delta
    lambda1 <- unlist(lambda0[names(theta.mle)])/abs(theta.mle)
	idx <- which(lambda>1e4)
	lambda[idx] <- 1e4 # lambda1[idx]
	
	f2 <- function( theta ) as.numeric( t(theta-theta.mle) %*% H %*% (theta-theta.mle) + lambda %*% abs(theta) )
	
	cat("\nPerforming LASSO estimation...\n")
	
	fit2 <- try( optim(theta.mle, f2, hessian=TRUE,..., 
					   control = list(maxit=30000, temp=2000, REPORT=500)), silent=TRUE) 
	if(class(fit2)=="try-error")
	return(list(mle=fail, sd.mle=NA, lasso=fail, sd.lasso=NA))
	
	theta.lasso <- fit2$par
	
	SIGMA1 <- try(sqrt(diag(solve(fit2$hessian))), silent=TRUE)
	
	if(class(SIGMA1)=="try-error")
	return(list(mle = theta.mle, sd.mle = NA, lasso = theta.lasso, sd.lasso = NA))
#	return(list(mle=fail, sd.mle=NA, lasso=fail, sd.lasso=NA))
	
	return(list(mle=theta.mle, sd.mle=SIGMA, lasso=theta.lasso, sd.lasso=SIGMA1,call=call, lambda0=lambda0))
}
