##################################################################
######  "Random number generators and 
######          related density functions for yuima packages"
######  (Last modified: March 9, 2010)
##################################################################


## Bilateral gamma

rbgamma <- function(x,delta.plus,gamma.plus,delta.minus,gamma.minus){
 if( delta.plus <= 0 )
   stop("delta.plus must be positive.")
 if( gamma.plus <= 0 )
   stop("gamma.plus must be positive.")
 if( delta.minus <= 0 )
   stop("delta.minus must be positive.")
 if( gamma.minus <= 0 )
   stop("gamma.minus must be positive.")
 X <- rgamma(x,delta.plus,gamma.plus) - rgamma(x,delta.minus,gamma.minus)
 return(X)
}


## (Multivariate) Normal gamma

rngamma <- function(x,lambda,alpha,beta,mu,Lambda){
  ## Error check
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda)){
    ## univariate case
    if(length(mu)!=1 || length(beta)!=1){
      stop("Error: wrong input dimension.")
    }
    tmp <- as.numeric(alpha^2 - beta^2)
    if( lambda <= 0 ){
      stop("lambda must be positive.")
    }
    if( alpha <= 0 ){
      stop("alpha must be positive.")
    }
    if( tmp <= 0 ){
      stop("alpha^2 - beta^2 must be positive value.")
    }
  
    tau <- rgamma(x,lambda,tmp/2)
    eta <- rnorm(x)
    ##  z <- mu + beta * tau * Lambda + sqrt(tau * Lambda) * eta
    z <- mu + beta * tau + sqrt(tau) * eta
    X <- z
    return(X)
    
  }else{ ## multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    if( sum(eigen(Lambda)$values <= 0 ) > 0){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) < 10^(-15)){
      stop("Determinant of Lambda must be one.")
    }
    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
  
    if( lambda <= 0 )
      stop("lambda must be positive.")
    if( alpha <= 0 )
      stop("alpha must be positive.")
    if( tmp <=0)
      stop("alpha^2 - t(beta) %*% Lambda %*% must be positive.")
    
    tau <- rgamma(x,lambda,tmp/2)
    eta <- rnorm(x*length(beta))
    sqrt.L <- svd(Lambda)
    sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)
    
    z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
    X <- z
    return(X)
  }
}


dngamma <- function(x,lam,al,be,mu){
    if( lam <= 0 )
	stop("lambda must be positive.")
    if( al <= 0 )
	stop("alpha must be positive.")
    if( al^2 - be^2 <= 0 )
	stop("alpha^2 - beta^2 must be positive.")
	
	dens <- exp(be*(x-mu))*((al^2 - be^2)^(lam))*besselK(al*abs(x-mu),lam-1/2)*abs(x-mu)^(lam-1/2)/(gamma(lam)*sqrt(pi)*((2*al)^(lam-1/2)))
	dens
}


## ## One-dim. normal gamma
## rngamma <- function(x,lambda,alpha,beta,mu){
##   tmp <- alpha^2 - beta^2
##   if( lambda <= 0 )
##     stop("lambda must be positive value.")
##   if( alpha <= 0 )
##     stop("alpha must be positive value.")
##   if( tmp <= 0 )
##     stop("alpha^2 - beta^2*mu must be positive value.")
## #  if( det(Lambda) != 1)
## ##  if( Lambda !=1)
## ##    stop("Determinant of Lambda must be one.")
  
##   tau <- rgamma(x,lambda,tmp/2)
##   eta <- rnorm(x)
## ##  z <- mu + beta * tau * Lambda + sqrt(tau * Lambda) * eta
##   z <- mu + beta * tau + sqrt(tau) * eta
##   X <- z
##   return(X)
## }

## rngamma <- function(x,lambda,alpha,beta,mu,Lambda){
##   tmp <- alpha^2 - t(beta) %*% Lambda %*% beta
##   if( lambda <= 0 )
##     stop("lambda must be positive value.")
##   if( alpha <= 0 )
##     stop("alpha must be positive value.")
##   if( tmp <= 0 )
##     stop("alpha^2 - beta^2*mu must be positive value.")
## #  if( det(Lambda) != 1)
##   if( Lambda !=1)
##     stop("Determinant of Lambda must be one.")
  
##   tau <- rgamma(x,lambda,tmp/2)
##   eta <- rnorm(x)
##   z <- mu + beta * tau * Lambda + sqrt(tau * Lambda) * eta
##   X <- z
##   return(X)
## }


## Inverse Gaussian

rIG <- function(x,delta,gamma){
  if( delta <= 0 )
    stop("delta must be positive.")
  if( gamma <= 0 )
    stop("gamma must be positive.")
  V <- rchisq(x,df=1);
  
  z1 <- ( delta/gamma + V/(2*gamma^2) ) - sqrt( V*delta/(gamma^3) + ( V/(2*gamma^2) )^2 )  
  U <- runif(x,min=0,max=1)
  idx <- which( U < (delta/(delta+gamma*z1)) )
  z2 <- (delta/gamma)^2 /z1[-idx]
  ret <- numeric(x)
  ret[idx] <- z1[idx]
  ret[-idx] <- z2
  
  return(ret)
}


dIG <- function(x,delta,gamma){
	if( delta <= 0 )
    stop("delta must be positive.")
	if( gamma <= 0 )
    stop("gamma must be positive.")

	dens <- delta*exp(delta*gamma)*(x^(-3/2))*exp(-((delta^2)/x+x*gamma^2)/2)/sqrt(2*pi)
	dens
}


## (Multivariate) Normal inverse Gaussian

rNIG <- function(x,alpha=1,beta=0,delta=1,mu=0,Lambda){
  ## Error check
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  
  if(missing(Lambda)){
    ## univariate case
    gamma <- sqrt(alpha^2 - beta^2)
    if(gamma <0){
      stop("alpha^2-beta^2 must be positive.")
    }
    
    if (gamma == 0) {
      V = rnorm(x)^2
      Z = delta * delta/V
      X = sqrt(Z) * rnorm(x)
    }else{ 
      Z <- rIG(x,delta,gamma)
      N <- rnorm(x,0,1)
      X <- mu + beta*Z + sqrt(Z)*N
    }
    return(X)
    
  }else{  ## multivariate case
    if( det(Lambda) < 10^(-15)){
      stop("Determinant of Lambda must be one.")
    }
    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
    if(tmp <0){
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative.")
    }
    gamma <- sqrt(tmp)
  }
  
  tau <- rIG(x,delta,gamma)
  eta <- rnorm(x*length(beta))
  sqrt.L <- svd(Lambda)
  sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)

  z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
  X <- z
  return(X)
}


dNIG <- function(x,al,be,de,mu){
	if( al < 0 )
	stop("alpha must be nonnegative.")
    if( al^2 - be^2 < 0 )
	stop("alpha^2 - beta^2 must be nonnegative.")
    if( de <= 0 )
	stop("delta must be positive.")

	dens <- al*de*exp(de*sqrt(al^{2}-be^{2})+be*(x-mu))*besselK(al*sqrt(de^{2}+(x-mu)^{2}),1)/(pi*sqrt(de^{2}+(x-mu)^{2}))
	dens
}


## ## One-dim. NIG
## rNIG <- function(x,alpha=1,beta=0,delta=1,mu=0){
##   gamma <- sqrt(alpha^2 - beta^2)
##   if (gamma == 0) {
##     V = rnorm(x)^2
##     Z = delta * delta/V
##     X = sqrt(Z) * rnorm(x)
##   }else{ 
##     Z <- rIG(x,delta,gamma)
##     N <- rnorm(x,0,1)
##     X <- mu + beta*Z + sqrt(Z)*N
##   }
##   return(X)
## }


## Non-Gaussian stable

rstable <- function(x,alpha,beta,sigma,gamma){
  a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
  b <- atan(beta*tan(alpha*pi/2))/alpha
  
  u <- runif(x,-pi/2,pi/2)
  v <- rexp(x,1)
  
  if(alpha!=1){
    s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
	X <- sigma * s + gamma * rep(1,x)
  }else{
    s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log(v*cos(u)/(beta*u + pi/2)))
	X <- sigma * s + (gamma + (2/pi)*beta*sigma*log(sigma)) * rep(1,x)
  }
  
  return(X)
}

## rstable <- function(x,alpha,beta,sigma, gamma, eps){
##   a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
##   b <- atan(beta*tan(alpha*pi/2))/alpha

##   u <- runif(x,-pi/2,pi/2)
##   v <- rexp(x,1)
  
##   if(alpha!=1){
##     s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
##   }else{
##     s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log(v*cos(u)/(beta*u + pi/2)))
##   }
  
##   X <- (eps^(1/alpha)) * sigma * s + gamma * eps * rep(1,x)
##   return(X)
## }


## Positive exponentially tempered stable (AR method)
## ## This must be re-coded later!! (temporarily, rather inefficient)
rpts <- function(x,al,a,b) {
    if( al <= 0 | al>= 1 )
    stop("al must lie in (0,1).")
    if( a <= 0 )
    stop("a must be positive value.")
    if( b <= 0 )
    stop("b must be positive value.")

	sig <- (-a*gamma(-al)*cos(pi*al/2))^(1/al)
	y <- c()
	i <- 1
	while (i <= x) {
		u <- runif(1)
		v <- rstable(1,al,1,sig,0)
		w <- exp(-b*v)
		if (u < w ){
			y <- append(y,v)
			i <- i+1
		}
	}
	return(y)
}


