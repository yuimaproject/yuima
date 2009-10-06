## random number generators for yuima packages

rIG <- function(x,delta,gamma){
  if( delta <= 0 )
    stop("delta must be positive value.")
  if( gamma <= 0 )
    stop("deltat must be positive value.")
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


rNIG <- function(x,alpha=1,beta=0,delta=1,mu=0){
  gamma <- sqrt(alpha^2 - beta^2)
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
}

rbgamma <- function(x,delta.plus,gamma.plus,delta.minus,gamma.minus){
 if( delta.plus <= 0 )
   stop("delta.plus must be positive value.")
 if( gamma.plus <= 0 )
   stop("gamma.plus must be positive value.")
 if( delta.minus <= 0 )
   stop("delta.minus must be positive value.")
 if( gamma.minus <= 0 )
   stop("gamma.minus must be positive value.")
 X <- rgamma(x,delta.plus,gamma.plus) - rgamma(x,delta.minus,gamma.minus)
 return(X)
}


## temporaly, desined only for univariate. 
rngamma <- function(x,lambda,alpha,beta,mu,Lambda){
  tmp <- alpha^2 - t(beta) %*% Lambda %*% beta
  if( lambda <= 0 )
    stop("lambda must be positive value.")
  if( alpha <= 0 )
    stop("alpha must be positive value.")
  if( tmp <= 0 )
    stop("alpha^2 - beta^2*mu must be positive value.")
#  if( det(Lambda) != 1)
  if( Lambda !=1)
    stop("Determinant of Lambda must be one.")

  tau <- rgamma(x,lambda,tmp/2)
  eta <- rnorm(x)
  z <- mu + beta * tau * Lambda + sqrt(tau * Lambda) * eta
  X <- z
  return(X)
}


rstable <- function(x,alpha,beta,sigma, gamma, eps){
  a <- (1 + (beta*tan(alpha*pi/2))^2)^(1/(2*alpha))
  b <- atan(beta*tan(alpha*pi/2))/alpha

  u <- runif(x,-pi/2,pi/2)
  v <- rexp(x,1)
  
  if(alpha!=1){
    s <- a * (sin(alpha*(u+b))/cos(u)^(1/alpha)) * (cos(u-alpha*(u+b))/v)^((1-alpha)/alpha)
  }else{
    s <- (2/pi) * ((pi/2 +beta*u)*tan(u) - beta * log(v*cos(u)/(beta*u + pi/2)))
  }
  
  X <- (eps^(1/alpha)) * sigma * s + gamma * eps * rep(1,x)
  return(X)
}
