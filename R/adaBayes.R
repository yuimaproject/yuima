## Bayes estimator to calculate theta1 and theta2
## arguments
## sdemod : sdeModel object
## h : sampling interval 
## prior1: prior distribution of theta1 (returned value is scalar)
## prior2: prior distribution of theta2 (returned value is scalar)
## domain1,domain2 : domain matrixes (or vectors) of theta1 and theta2 (length(theta1) x 2 matrix ,length(theta2) x 2 matrix)
## ex.) theta1 <- c("theta11","theta12") , theta2 = c("theta21","theta22","theta23")
##     theta1.lower <- c(0,0) , theta1.upper <- c(1,2)
##     theta2.lower <- c(0.5,0.2,0.7) , theta2.upper <- c(10,5,15)
##     domain1 <- matrix(c(theta1.lower , theta1.upper ),2,2)
##     domain2 <- matrix(c(theta2.lower , theta2.upper ),3,2)
## theta1.ans , theta2.ans : offset vector (or scalar) to calculate
setGeneric("adaBayes",
           function(yuima,h,prior2,prior1,domain2,domain1,theta2.offset=0,theta1.offset=0)
           standardGeneric("adaBayes")
           )
setMethod("adaBayes",signature("yuima"),function(yuima,h,prior2,prior1,domain2 ,domain1 ,theta2.offset=0,theta1.offset=0){
  ## Hn function
  Hn <- function(yuima,h,theta1,theta2){
    d.size <- yuima@model@equation.number
    n.size <- dim(as.matrix(onezoo(yuima)))[1]
    tmp <- ((n.size*d.size/2) * log( (2*pi*h) ) + ql(yuima,theta2=theta2,theta1=theta1,h=h))
    return(tmp)
  }
  
  ## rHn function
  rHn <- function(yuima,h,theta1,theta2,theta2.offset,theta1.offset){
    d.size <- yuima@model@equation.number
    n.size <- dim(as.matrix(onezoo(yuima)))[1]
    tmp <- ((n.size*d.size/2) * log( (2*pi*h) ) + rql(yuima,theta2=theta2,theta1=theta1,theta2.offset,theta1.offset,h=h))
    return(tmp)
  }
  
  ## Bayes estimator for theta1
  bayes.theta1.tilda <- function(yuima,prior1=prior1){
    rHn.tmp <- function(yuima,h,theta1,theta2,theta2.offset,theta1.offset){
      tmp <- numeric( length(theta1) )
      for(i in 1:length(theta1)){
        rHn1 <- rHn( yuima , h , theta1[i] , theta2 ,theta2.offset,theta1.offset)
        tmp[i] <- exp( rHn1 )
      }
      if( is.infinite(sum(tmp)) != 0 ){
        stop("Hn is too big.\n Please check theta2.offset and theta1.offset")
      }
      return(tmp)
    }
    ## if theta1 is a vector , prior1 returns vector compliant with each theta1 value
    prior1.tmp <- function(theta1){
      tmp <- numeric(length(theta1))
      for( i in 1:length(theta1) ){
        tmp[i] <- prior1(theta1[i],domain1)
      }
      return(tmp)
    }
    ## numerator function 
    Numerator.tmp <- function(theta1){
      tmp <- theta1 *rHn.tmp(yuima=yuima,h=h,theta1,theta2=1,theta2.offset=theta2.offset,theta1.offset=theta1.offset) * prior1.tmp(theta1)
      return(tmp)
    }
    ## denominator function
    Denominator.tmp <- function(theta1){
      tmp <- rHn.tmp(yuima=yuima,h=h,theta1,theta2=1,theta2.offset=theta2.offset,theta1.offset=theta1.offset) * prior1.tmp(theta1)
      return(tmp)
    }
    numerator <- integrate(Numerator.tmp,domain1[1,1],domain1[1,2])
    denominator <- integrate(Denominator.tmp,domain1[1,1],domain1[1,2])
    if(numerator$value == 0 || denominator$value == 0){
      stop("Quasi-likelihood is too small.\n  Please check domain1 , theta2.offset and theta1.offset")
    }
    return( numerator$value/denominator$value )
  }
  
  ##Bayes estimator for theta2
  bayes.theta2.tilda <- function(yuima,prior2=prior2){
    Hn.tmp <- function(yuima,h,theta1,theta2){
      tmp <- numeric(length(theta2))
      for(i in 1:length(theta2)){
        Hn1 <- Hn( yuima , h , theta1 , theta2[i] )
        tmp[i] <- exp( Hn1 -Hn2 )      
      }
      if( is.infinite(sum(tmp)) != 0 ){
        stop("Hn is too big.\n  Please check theta2.offset and theta1.offset")
      }
      return(tmp)
    }
    ## if theta2 is a vector , prior1 returns vector compliant with each theta2 value
    prior2.tmp <- function(theta2){
      tmp <- numeric(length(theta2))
      for( i in 1:length(theta2) ){
        tmp[i] <- prior2(theta2[i],domain2)
      }
      return(tmp)
    }
    ## numerator function 
    Numerator.tmp <- function(theta2){
      tmp <- theta2 * Hn.tmp(yuima=yuima,h=h,theta1=theta1.tilda,theta2) * prior2.tmp(theta2)
      return(tmp)
    }
    ## denominator function
    Denominator.tmp <- function(theta2){
      tmp <- Hn.tmp(yuima=yuima,h=h,theta1=theta1.tilda,theta2) * prior2.tmp(theta2)
      return(tmp)
    }
    numerator <- integrate(Numerator.tmp,domain2[1,1],domain2[1,2])
    denominator <- integrate(Denominator.tmp,domain2[1,1],domain2[1,2])
    if(numerator$value == 0 || denominator$value == 0){
      stop("Quasi-likelihood is too small.\n  Please check domain2 , theta2.offset and theta1.offset")
    }
    return( numerator$value/denominator$value )
  }
  ## Start error Check ##
  ## check yuima
  if( missing(yuima) ){
    stop("yuima is missing.")
  }
  ## check number of theta
  if( length(yuima@model@diffusion) != 1 || length(yuima@model@drift) != 1){
    stop("This version do not support multi-dimention parameters.")
  }else{
    if( is.vector(domain2) ) domain2 <- t(as.matrix(domain2))
    if( is.vector(domain1) ) domain1 <- t(as.matrix(domain1))
  }
  ## check domain
  if( !is.matrix(domain2) || !is.matrix(domain1) ){
    stop("domain2 or domain1 is not a matrix.")
  }else if( length(yuima@model@diffusion) != nrow(domain1) || length(yuima@model@drift) != nrow(domain2)){
    stop("dimention of domain1 or domain2 is different from length with respect to each theta.")
  }else if( ncol(domain1) != 2 || ncol(domain2) != 2){
    stop("ncol(domain1) or ncol(domain2) is irrelevance.")
  }
  ## check X
      ##if(is.matrix(yuima@data@original.data) == FALSE){##original.data
      ##  stop("data is not a matrix. Please check again.")
      ##}else if
  if( dim(as.matrix(onezoo(yuima)))[1] < 2 ){
    stop("length of data is too short.")
  }
  ## check h
  if( (!is.double(h) && !is.integer(h)) || length(h) != 1){
    stop("'h' is not correct type.")
  }else if( h <= 0 ){
    stop("'h' is negative or zero. Please set 'h' to be a positive value.")
  }
  ## check theta offset
  if( !is.null(theta1.offset) && !is.null(theta2.offset) ){
    if(!is.vector(theta2.offset) || !is.vector(theta1.offset)){
      stop("theta2.offset or theta1.offset is not a vector.")
    }
    if(length(yuima@model@drift) != length(theta2.offset) || length(yuima@model@diffusion) != length(theta1.offset)){
      stop("length of theta2.offset or theta1.offset is different from length with respect to each theta")
    }
    Hn2 <- Hn( yuima , h , theta1 = theta1.offset , theta2 = theta2.offset  )
  }else{
    Hn2 <- 0
  }
  ## Finish error check ##
  
  ## Culculate theta1.tilda and theta2.tilda
  theta1.tilda <- bayes.theta1.tilda(yuima,prior1=prior1)
  theta2.tilda <- bayes.theta2.tilda(yuima,prior2=prior2)
  return(c(theta2.tilda,theta1.tilda))
}
)
