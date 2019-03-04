##################################################################
######  "Random number generators and 
######          related density functions for yuima packages"
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

# "dbgamma" by YU.
dbgamma<-function(x,delta.plus,gamma.plus,delta.minus,gamma.minus){
  ## Error check
  if(length(delta.plus)!=1||length(gamma.plus)!=1||length(delta.minus)!=1||length(gamma.minus)!=1)
  stop("All of the parameters are numeric.")
  if(delta.plus<=0||gamma.plus<=0||delta.minus<=0||gamma.minus<=0)
  stop("All of the parameters are positive.")
  
  leng <- length(x)
  dens <- numeric(leng)
  
  for(i in 1:leng){
  if(x[i]>=0){
  ## On the positive line
  funcp<-function(x,y){y^{delta.minus-1}*(x+y/(gamma.plus+gamma.minus))^{delta.plus-1}*exp(-y)}
  intp<-function(x){integrate(funcp,lower=0,upper=Inf,x=x)$value}
  intvecp<-function(x)sapply(x,intp)
  dens[i]<-gamma.plus^delta.plus*gamma.minus^delta.minus/((gamma.plus+gamma.minus)^delta.minus*gamma(delta.plus)*gamma(delta.minus))*exp(-gamma.plus*x[i])*intvecp(x[i])
  }else{
  ## On the negative line
  funcm<-function(x,y){y^{delta.plus-1}*(-x+y/(gamma.plus+gamma.minus))^{delta.minus-1}*exp(-y)}
  intm<-function(x){integrate(funcm,lower=0,upper=Inf,x=x)$value}
  intvecm<-function(x)sapply(x,intm)
  dens[i]<-gamma.plus^delta.plus*gamma.minus^delta.minus/((gamma.plus+gamma.minus)^delta.plus*gamma(delta.minus)*gamma(delta.plus))*exp(gamma.minus*x[i])*intvecm(x[i])
  }
  }
  dens
}

## Generalized inverse Gaussian

rGIG<-function(x,lambda,delta,gamma){
    if((delta < 0)||(gamma<0)){
        stop("delta and gamma must be nonnegative.")
    }
    if((delta == 0) && (lambda <= 0)){
        stop("delta must be positive when lambda is nonpositive.")
    }
    if((gamma == 0) && (lambda >= 0)){
        stop("gamma must be positive when lambda is nonnegative.")
    }
    if(lambda >= 0){
        if(delta == 0){ ## gamma case
            X<-rgamma(x, lambda, 1/(2*gamma^2))
        }
        if(gamma == 0){ ## inverse gamma case
            X<-delta^2/(2*rgamma(x,-lambda,1))
        }
        if(delta != 0){
            X<-.C("rGIG", as.integer(x), as.double(lambda), as.double(delta^2), as.double(gamma^2),
            rn = double(length = x))$rn}
    }
    else{
        Y<-.C("rGIG", as.integer(x), as.double(-lambda), as.double(gamma^2), as.double(delta^2),
        rn = double(length = x))$rn
        X<-1/Y
    } 
    return(X)
}

dGIG<-function(x,lambda,delta,gamma){
    if(x <= 0)
    stop("x must be positive")
    if((delta < 0)||(gamma<0)){
        stop("delta and gamma must be nonnegative.")
    }
    if((delta == 0) && (lambda <= 0)){
        stop("delta must be positive when lambda is nonpositive.")
    }
    if((gamma == 0) && (lambda >= 0)){
        stop("gamma must be positive when lambda is nonnegative.")
    }
    if(delta == 0){ ## gamma case
        dens <- dgamma(x,lambda,1/(2*gamma^2))
    }
    if(gamma == 0){ ## inverse gamma case
        dens <- x^(lambda-1)*exp(-1/2*delta/x)/(gamma(-lambda)*2^(-lambda)*delta^(2*lambda))
    }else{
        dens <- 1/2*(gamma/delta)^lambda*1/besselK(gamma*delta,lambda)*x^(lambda-1)*exp(-1/2*(delta^2/x+gamma^2*x))
    }
    return(dens)
}

## Generalized hyperbolic

rGH<-function(x,lambda,alpha,beta,delta,mu,Lambda){
    if(length(mu)!=length(beta)){
        stop("Error: wrong input dimension.")
    }
    
    if(missing(Lambda))
    Lambda <- NA
    
    if( alpha < 0 ){
        stop("alpha must be nonnegative.")
    }
    
    ## variance gamma case
    if(delta == 0){
        X<-rvgamma(x,lambda,alpha,beta,mu,lambda)
    }
    
    ## univariate case
    if(any(is.na(Lambda))){
        if(length(mu)!=1 || length(beta)!=1){
            stop("Error: wrong input dimension.")
        }
        
        tmp <- as.numeric(alpha^2 - beta^2)
        if( tmp < 0 ){
            stop("alpha^2 - beta^2 must be nonnegative value.")
        }
        
        if(((alpha == 0)||(tmp == 0)) && (lambda >= 0)){
            stop("alpha and alpha^2 - beta^2 must be positive when lambda is nonnegative.")
        }
        
        Z<-rGIG(x,lambda,delta,sqrt(tmp))
        N<-rnorm(x,0,1)
        X<-mu + beta*Z + sqrt(Z)*N
    }
    else{## multivariate case
        if( nrow(Lambda)!=ncol(Lambda)){
            stop("Lambda must be a square matrix.")
        }
        
        if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
            stop("Lambda must be a symmetric matrix")
        }
        
        if( nrow(Lambda)!=length(beta)){
            stop("Dimension of Lambda and beta must be equal.")
        }
        
        if( min(eigen(Lambda)$value) <= 10^(-15) ){
            stop("Lambda must be positive definite.")
        }
        
        if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
            stop("The determinant of Lambda must be 1.")
        }
        
        tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
        
        if(tmp < 0){
            stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative")
        }
        
        if(((alpha == 0)||(tmp == 0))&&(lambda >=0)){
            stop("alpha and alpha^2 - t(beta) %*% Lambda %*% beta must be positive when lambda is nonnegative.")
        }
        
        Z<-rGIG(x,lambda,delta,sqrt(tmp))
        N<-rnorm(x*length(beta),0,1)
        sqrt.L <- svd(Lambda)
        sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)
        X <- mu + t(matrix(rep(Z,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(Z),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(N,x,length(beta))))
    }
    return(X)
}


dGH<-function(x,lambda,alpha,beta,delta,mu,Lambda){
    if(length(mu)!=length(beta)){
        stop("Error: wrong input dimension.")
    }
    
    if(missing(Lambda))
    Lambda <- NA
    
    if( alpha < 0 ){
        stop("alpha must be nonnegative.")
    }
    
    ## variance gamma case
    if(delta == 0){
        X<-dvgamma(x,lambda,alpha,beta,mu,Lambda)
    }
    
    ## univariate case
    if(any(is.na(Lambda))){
        if(length(mu)!=1 || length(beta)!=1){
            stop("Error: wrong input dimension.")
        }
        
        tmp <- as.numeric(alpha^2 - beta^2)
        if( tmp < 0 ){
            stop("alpha^2 - beta^2 must be nonnegative value.")
        }
        
        if(((alpha == 0)||(tmp == 0)) && (lambda >= 0)){
            stop("alpha and alpha^2 - beta^2 must be positive when lambda is nonnegative.")
        }
        
        nu<--2*lambda
        
        if(alpha == 0){## gamma = 0 (in gig), scaled t-distribution
            dens<-gamma((nu+1)/2)*(1+((x-mu)/delta)^2)^(-(nu+1)/2)/(delta*sqrt(pi)*gamma(nu/2))
        }else if(tmp == 0){## skewed t-distribution
            dens<-delta^nu*exp(beta*(x-mu))*gamma((nu+1)/2)*besselK(abs(beta)*sqrt(delta^2+(x-mu)^2),(nu+1)/2)/(2^((nu-1)/2)*sqrt(pi)*gamma(nu/2)*(sqrt(delta^2+(x-mu)^2)/abs(beta))^((nu+1)/2))
        }else{
            dens<-tmp^(lambda/2)*sqrt(delta^2+(x-mu)^2)^(lambda-1/2)*besselK(alpha*sqrt(delta^2+(x-mu)^2),lambda-1/2)*exp(beta*(x-mu))/(sqrt(2*pi)*alpha^(lambda-1/2)*delta^lambda*besselK(delta*sqrt(tmp),lambda))
        }
    }
    else{## multivariate case
        if( nrow(Lambda)!=ncol(Lambda)){
            stop("Lambda must be a square matrix.")
        }
        
        if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
            stop("Lambda must be a symmetric matrix")
        }
        
        if( nrow(Lambda)!=length(beta)){
            stop("Dimension of Lambda and beta must be equal.")
        }
        
        if( min(eigen(Lambda)$value) <= 10^(-15) ){
            stop("Lambda must be positive definite.")
        }
        
        if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
            stop("The determinant of Lambda must be 1.")
        }
        
        tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
        
        if(tmp < 0){
            stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative")
        }
        
        if(((alpha == 0)||(tmp == 0))&&(lambda >=0)){
            stop("alpha and alpha^2 - t(beta) %*% Lambda %*% beta must be positive when lambda is nonnegative.")
        }
        
        Lambdainv<-solve(Lambda)
        nu<--2*lambda
        d<-nrow(Lambda)
        
        if(alpha == 0){ ## gamma = 0 (in gig) multivariate scaled t-distribution
            dens<-gamma((nu+d)/2)*(1+t(x-mu)%*%Lambdainv%*%(x-mu)/delta^2)^(-(nu+d)/2)/(pi^(d/2)*gamma(nu/2)*delta^d)
        }else if(tmp == 0){ ## multivariate skewed t-distribution
            dens<-delta^nu*exp(t(beta)%*%(x-mu))*besselK(alpha*sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)),-(nu+d)/2)/((pi)^(d/2)*2^((nu+d)/2-1)*(sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu))/alpha)^((nu+d)/2)*gamma(nu/2))
        }else{
            dens<-exp(t(beta)%*%(x-mu))*(sqrt(tmp)/delta)^lambda*besselK(alpha*sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)),-(nu+d)/2)/((2*pi)^(d/2)*besselK(delta*sqrt(tmp),lambda)*(sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu))/alpha)^(d/2-lambda))
        }
    }
    return(dens)
}


## (Multivariate) Variance gamma

rvgamma <- function(x,lambda,alpha,beta,mu,Lambda){
  ## Error check
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
   Lambda <- NA
  
  if(any(is.na(Lambda))){
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
    if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
      stop("Lambda must be a symmetric matrix")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }

    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
  
    if( lambda <= 0 )
      stop("lambda must be positive.")
    if( alpha <= 0 )
      stop("alpha must be positive.")
    if( tmp <=0)
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be positive.")
    
    tau <- rgamma(x,lambda,tmp/2)
    eta <- rnorm(x*length(beta))
    sqrt.L <- svd(Lambda)
    sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)
    
    z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
    X <- z
    return(X)
  }
}


dvgamma <- function(x,lambda,alpha,beta,mu,Lambda){
  ## Error check
  if(length(lambda)!=1||length(alpha)!=1)
    stop("alpha and lambda must be positive reals.")
  if( lambda <= 0 )
    stop("lambda must be positive.")
  if( alpha <= 0 )
    stop("alpha must be positive.")
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
    Lambda <- NA
  if(any(is.na(Lambda))){
    ## univariate case
    if(length(mu)!=1 || length(beta)!=1){
      stop("Error: wrong input dimension.")
    }
    if( alpha^2 - beta^2 <= 0 )
      stop("alpha^2 - beta^2 must be positive.")
    
    dens <- exp(beta*(x-mu))*((alpha^2 - beta^2)^(lambda))*besselK(alpha*abs(x-mu),lambda-1/2)*abs(x-mu)^(lambda-1/2)/(gamma(lambda)*sqrt(pi)*((2*alpha)^(lambda-1/2)))
    dens}
  else{ ## multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
      stop("Lambda must be a symmetric matrix")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    
    if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }

    tmp <- as.numeric(alpha^2 - t(beta) %*% Lambda %*% beta)
    if( tmp <=0)
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be positive.")
    Lambdainv<-solve(Lambda)
    dens<- exp(t(beta)%*%(x-mu))*(alpha^2-t(beta)%*%Lambda%*%beta)^(lambda)*besselK(alpha*sqrt(t(x-mu)%*%Lambdainv%*%(x-mu)),lambda-nrow(Lambda)/2)*sqrt(t(x-mu)%*%Lambdainv%*%(x-mu))^{lambda-nrow(Lambda)/2}/(gamma(lambda)*pi^{nrow(Lambda)/2}*2^{nrow(Lambda)/2+lambda-1}*alpha^{lambda-nrow(Lambda)/2})
    dens
  }
}


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
  if(x <= 0)
    stop("x must be positive")
	if( delta <= 0 )
    stop("delta must be positive.")
	if( gamma <= 0 )
    stop("gamma must be positive.")

	dens <- delta*exp(delta*gamma)*(x^(-3/2))*exp(-((delta^2)/x+x*gamma^2)/2)/sqrt(2*pi)
	dens
}


## (Multivariate) Normal inverse Gaussian

rNIG <- function(x,alpha,beta,delta,mu,Lambda){
  ## Error check
  #print(match.call())
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(length(alpha)>1||length(delta)>1)
    stop("alpha and delta must be positive reals.")
  if( alpha < 0 )
    stop("alpha must be nonnegative.")
  if( delta <= 0 )
    stop("delta must be positive.")
  
  if(missing(Lambda))
   Lambda <- NA

  if(any(is.na(Lambda)) & length(Lambda)==1){
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
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
      stop("Lambda must be a symmetric matrix")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
	if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
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
  #  print(str(X))
  return(X)
}


dNIG <- function(x,alpha,beta,delta,mu,Lambda){
  ## Error check
  #print(match.call())
  if(length(alpha)>1||length(delta)>1)
    stop("alpha and delta must be positive reals.")
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if( alpha < 0 )
    stop("alpha must be nonnegative.")
  if( delta <= 0 )
    stop("delta must be positive.")
  if(missing(Lambda))
    Lambda <- NA
  if(any(is.na(Lambda))){
    #univraiate case
    if(length(beta)>1||length(mu)>1)
      stop("beta and mu must be numeric")
    if( alpha^2 - beta^2 < 0 )
      stop("alpha^2 - beta^2 must be nonnegative.")
    dens <- alpha*delta*exp(delta*sqrt(alpha^{2}-beta^{2})+beta*(x-mu))*besselK(alpha*sqrt(delta^{2}+(x-mu)^{2}),1)/(pi*sqrt(delta^{2}+(x-mu)^{2}))
    dens
  }else{
    #multivariate case
    if( nrow(Lambda)!=ncol(Lambda)){
      stop("Lambda must be a square matrix.")
    }
    if(sum((Lambda-t(Lambda))*(Lambda-t(Lambda)))!=0){
      stop("Lambda must be a symmetric matrix")
    }
    if( nrow(Lambda)!=length(beta)){
      stop("Dimension of Lambda and beta must be equal.")
    }
    if(nrow(Lambda)!=length(mu)){
      stop("Dimension of Lambda and mu must be equal.")
    }

	if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
    	stop("The determinant of Lambda must be 1.")
    }
    
    tmp<-alpha-t(beta)%*%Lambda%*%beta
    if(tmp <0){
      stop("alpha^2 - t(beta) %*% Lambda %*% beta must be nonnegative.")
    }
    Lambdainv<-solve(Lambda)
    dens<- alpha^{(nrow(Lambda)+1)/2}*delta*exp(delta*sqrt(tmp)+t(beta)%*%(x-mu))*besselK(alpha*sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)),nrow(Lambda))/(2^{(nrow(Lambda)-1)/2}*pi^{(nrow(Lambda)+1)/2}*(sqrt(delta^2+t(x-mu)%*%Lambdainv%*%(x-mu)))^{(nrow(Lambda)+1)/2})
    dens
  }	
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


## Univariate non-Gaussian stable: corrected Weron's algorithm incorporated
## (20160114, HM) "dstable" still unavailable in YUIMA... incorporate "stabledist" package?

rstable <- function(x,alpha,beta,sigma,gamma){
	if( alpha <= 0 || alpha > 2)
	stop("The index alpha must take values in (0,2].")
	if( beta < -1 || beta > 1)
	stop("The skeweness beta must take values in [-1,1].")
	if( sigma <= 0 )
	stop("The scale sigma must be positive.")
	
  a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
  b <- atan(beta*tan(alpha*pi/2))/alpha
  
  u <- runif(x,-pi/2,pi/2)
  v <- rexp(x,1)
  
  if(alpha!=1){
    s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
	X <- sigma * s + gamma * rep(1,x)
  }else{
    s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log((pi/2)*v*cos(u)/(beta*u + pi/2)))
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
rpts<-function(x,alpha,a,b){
  if( alpha <= 0 | alpha>= 1 )
    stop("alpha must lie in (0,1).")
  if( a <= 0 )
    stop("a must be positive value.")
  if( b <= 0 )
    stop("b must be positive value.")
  ar<-exp(a*gamma(-alpha)*b^(alpha))
  if(ar <= 0.1)
    stop("Acceptance rate is too small.")
  else
    .C("rpts",as.integer(x),as.double(alpha),as.double(a),as.double(b),rn=double(length=x))$rn}


## Normal tempered stable 
rnts<-function(x,alpha,a,b,beta,mu,Lambda){
  ## Error check
  if(length(mu)!=length(beta)){
    stop("Error: wrong input dimension.")
  }
  if(missing(Lambda))
    Lambda <- NA
  if( alpha <= 0 || alpha > 2 ){
    stop("The index alpha must take values in (0,2].")
  }
  if( a <= 0 )
    stop("a must be positive value.")
  if( b <= 0 )
    stop("b must be positive value.")
  
  if(any(is.na(Lambda))){
    ## univariate case
    if(length(mu)!=1 || length(beta)!=1){
      stop("Error: wrong input dimension.")
    }
    
    tau <- rpts(x,alpha,a,b)
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
    if(nrow(Lambda)!=length(mu)){
      stop("Dimension of Lambda and mu must be equal.")
    }
    if( min(eigen(Lambda)$value) <= 10^(-15) ){
      stop("Lambda must be positive definite.")
    }
    if( det(Lambda) > 1+10^(-15) || det(Lambda) < 1-10^(-15) ){
      stop("The determinant of Lambda must be 1.")
    }
    tau<-rpts(x,alpha,a,b)
    eta <- rnorm(x*length(beta))
    sqrt.L <- svd(Lambda)
    sqrt.L <- sqrt.L$u %*% diag(sqrt(sqrt.L$d)) %*% t(sqrt.L$v)
    
    z <- mu + t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
    z<-mu+ t(matrix(rep(tau,length(beta)),x,length(beta))) * matrix(rep(Lambda %*% beta,x),length(beta),x)+t(matrix(rep(sqrt(tau),length(beta)),x,length(beta))) * (sqrt.L %*% t(matrix(eta,x,length(beta))))
    X <- z
    return(X)
  }
}
