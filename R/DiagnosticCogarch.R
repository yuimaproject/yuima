StartValCog<-function(yuima.cogarch, param=list(), mu=1, rho=3){
  
  
  if(missing(yuima.cogarch))
    yuima.stop("yuima.cogarch or yuima object is missing.")
  
  if(length(param)==0)
    yuima.stop("missing values parameters")
    
  if(!is.COGARCH(yuima.cogarch)){
    yuima.warn("The model does not belong to the class yuima.cogarch")
  }
  
  
  if(is(yuima.cogarch,"yuima")){
    model<-yuima.cogarch@model
  }else{
    if(is(yuima.cogarch,"yuima.cogarch")){
      model<-yuima.cogarch
    }
  }
  
  info <- model@info
  numb.ar <- info@q
  ar.name <- paste(info@ar.par,c(numb.ar:1),sep="")
  numb.ma <- info@p
  ma.name <- paste(info@ma.par,c(1:numb.ma),sep="")
  loc.par <- info@loc.par
  
  nm <- c(names(param))
  param<-as.numeric(param)
  names(param)<-nm
  
  
  xinit.name0 <- model@parameter@xinit
  idx <- na.omit(match(c(loc.par, ma.name), xinit.name0))
  xinit.name <- xinit.name0[-idx]
  fullcoeff <- c(ar.name, ma.name, loc.par,xinit.name)
  
  oo <- match(nm, fullcoeff)
  
  if(any(is.na(oo)))
    yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima.cogarch model")
 
  acoeff <- param[ma.name]
  b <- param[ar.name]
  cost<- param[loc.par]
  res<-StationaryMoments(cost,b,acoeff,mu,rho)
  return(res)

}

StationaryMoments<-function(cost,b,acoeff,mu=1,rho=3){
  # We obtain stationary mean of State process
  # stationary mean of the variance 
  # E\left(Y\right)=-a_{0}m_{2}\left(A+m_{2}ea'\right)^{-1}e
  q<-length(b) 
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  p<-length(acoeff)
  a[1:p,1] <- acoeff
  B_tilde <- MatrixA(b[c(q:1)])+mu*e%*%t(a)
  
  if(q>1){
    invB<-rbind(c(-B_tilde[q,-1],1)/B_tilde[q,1],cbind(diag(q-1),matrix(0,q-1,1)))
  }else{invB<-1/B_tilde}
  
  ExpStatVar <- -cost*mu*invB%*%e
  ExpVar <- cost+t(a)%*%ExpStatVar
  res <- list(ExpVar=ExpVar, ExpStatVar=ExpStatVar)
  return(res)
}

