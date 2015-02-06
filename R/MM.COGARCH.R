# In this function we consider different kind estimation procedures for COGARCH(P,Q)
# Model. 
is.COGARCH <- function(obj){
  if(is(obj,"yuima"))
    return(is(obj@model, "yuima.cogarch"))
  if(is(obj,"yuima.cogarch"))
    return(is(obj, "yuima.cogarch"))
  return(FALSE)
}

# The estimation procedure for cogarch(p,q) implemented in this code are based on the 
# Chadraa phd's thesis
MM.COGARCH<-function(yuima, data = NULL, start, method="BFGS", fixed = list(), 
                           lower, upper){
  print <- FALSE
  call <- match.call()
# First step we check if all inputs are inserted correctly.
  if( missing(yuima))
    yuima.stop("yuima object is missing.")
  
  if( missing(start) ) 
    yuima.stop("Starting values for the parameters are missing.")
  
  if( !is.COGARCH(yuima) )
    yuima.stop("The model is not an object of class yuima.")

  if( !is(yuima,"yuima") && missing(data) )
    yuima.stop("data are missing.")
  
  
  if(is(yuima,"yuima")){
    model<-yuima@model
    if(is.null(data)){
      observ<-yuima@data
    }else{observ<-data}
  }else{
    if(is(yuima,"yuima.cogarch")){
      model<-yuima
      if(is.null(data)){
        yuima.stop("Missing data")
      }
      observ<-data
    }
  }
  
  if(!is(observ,"yuima.data")){
    yuima.stop("Data must be an object of class yuima.data-class")  
  } 
  
  if( !missing(upper)){
    yuima.warn("The upper constraints will be implemented as soon as possible")
    upper <- list()
  }

  if( !missing(lower)){
    yuima.warn("The lower constraints will be implemented as soon as possible")
    lower <- list()
  }

  if( !missing(fixed)){
    yuima.warn("The fixed constraints will be implemented as soon as possible")
    fixed <- list()
  }

  # We identify the model parameters
  info <- model@info
  
  numb.ar <- info@q
  ar.name <- paste(info@ar.par,c(numb.ar:1),sep="")
  
  numb.ma <- info@p
  ma.name <- paste(info@ma.par,c(1:numb.ma),sep="")
  
  loc.par <- info@loc.par
  xinit.name <- paste(info@Latent.var, c(1:numb.ar), sep = "")

  meas.par <- model@parameter@measure
  
  fixed.name <- names(fixed)
  if(info@q==1){
  #  nm <- c(names(start), "EL1", "phi1","phi2")
    nm <- c(names(start), "EL1")
  }else{
  #  nm <- c(names(start), "EL1", "M2Lev","M4Lev")
    nm <- c(names(start), "EL1")
  }
  # We identify the index of parameters
  if(length(meas.par)!=0){
    if(info@q==1){  
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, measure.par, "EL1", "phi1","phi2")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, meas.par, "EL1")
    }else{
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, measure.par, "EL1", "M2Lev","M4Lev")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name, meas.par, "EL1")
    }
  }else{
    if(info@q==1){
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1", "phi1","phi2")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1")
    }else{
  #    fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1", "M2Lev","M4Lev")
      fullcoeff <- c(ar.name, ma.name, loc.par, xinit.name,  "EL1")
    }
  }
  oo <- match(nm, fullcoeff)

  if(any(is.na(oo)))
    yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
  if(info@q==1){
  #  start <- c(start, EL1 = 1, phi1=-1, phi2=-1)
    start <- c(start, EL1 = 1)
  }else{
  #  start <- c(start, EL1 = 1, M2Lev=1, M4Lev=2)
    start <- c(start, EL1 = 1)
  }
  start <- start[order(oo)]
  
  nm <- names(start)

  ar.idx <- match(ar.name, fullcoeff)
  ma.idx <- match(ma.name, fullcoeff)
  loc.idx <- match(loc.par, fullcoeff)
  meas.idx <- match(meas.par, fullcoeff)
  fixed.idx <- match(fixed.name, fullcoeff)
  EL1.idx <- match("EL1",fullcoeff)  # We decide to pass EL1 as parameter !!!
#   if(info@q==1){
#     phi1.idx <- match("phi1",fullcoeff) # We decide to pass phi1 as parameter !!! 
#     phi2.idx <- match("phi2",fullcoeff) # We decide to pass phi1 as parameter !!!
#   }else{
#     M2Lev.idx <- match("M2Lev",fullcoeff)
#     M4Lev.idx <- match("M4Lev",fullcoeff)
#   }
  # We build the environment used as input in the 
  # objective function.

  env <- new.env()
  n <- length(observ)[1]
  
  # Data
  assign("Data",  as.matrix(onezoo(observ)[,1]), envir=env)
  assign("deltaData",  frequency(onezoo(observ)[,1]), envir=env)
  assign("time.obs",length(env$Data),envir=env)
  
  # Order
  assign("p", info@p, envir=env)
  assign("q", info@q, envir=env)

  # Idx
  assign("ar.idx", ar.idx, envir=env)
  assign("ma.idx", ma.idx, envir=env)
  assign("loc.idx", loc.idx, envir=env)
  assign("meas.idx", meas.idx, envir=env)
  assign("EL1.idx", EL1.idx, envir=env)
#   if(info@q==1){ 
#     assign("phi1.idx", phi1.idx, envir=env)
#     assign("phi2.idx", phi1.idx, envir=env)
#   }else{
#     assign("M2Lev.idx", M2Lev.idx, envir=env) 
#     assign("M4Lev.idx", M4Lev.idx, envir=env)
#   }
    
objectiveFun <- function(p) {
  mycoef <- as.list(p)
  #		 names(mycoef) <- nm
  
    if(length(c(fixed.idx, meas.idx))>0){ ## SMI 2/9/14
      names(mycoef) <- nm[-c(fixed.idx,meas.idx)] ## SMI 2/9/14
    }else{
      names(mycoef) <- nm
   } 
  ErrTerm(yuima=yuima, param=mycoef, print=print, env)
}
 out<- optim(start, objectiveFun, method=method)
# , method = "L-BFGS-B",
#              lower=c(0.0001,0.00001,0.00001,0, 0.9,-100,-100),
#              upper=c(1,1,1,1,1.1,0,0))

  return(out)
  
}

ErrTerm <- function(yuima, param, print, env){
 # sample moments
  typeacf <- "correlation"
  param <- as.numeric(param)
  G_i <- diff(env$Data)
  r <- 1/env$deltaData 
  mu_G2 <- mean(G_i^2)
  var_G2 <- mean(G_i^4) - mu_G2^2
  d <- floor(sqrt(length(G_i)))
  CovQuad <- log(abs(acf(G_i^2,plot=FALSE,lag.max=min(d,10),type=typeacf)$acf[-1]))
  h <- seq(1, d, by = 1)*r
  cost <- env$loc.idx
  b <- env$ar.idx
  a <- env$ma.idx
  meanL1 <- param[env$EL1.idx]
  #meanL1 <- 1
#   if(env$q == 1){
# 
#   beta <- param[cost]*param[b]
#   eta <- param[b]
#   phi <- param[a]
#    
# #   phi1 <- param[env$phi1.idx]
# #   phi2 <- param[env$phi2.idx]
#   #theo_mu_G2 <- meanL1*r*beta/abs(phi1)
#   phi1 <- meanL1*r*beta/mu_G2
#   
#   termA <- (6*mu_G2/r*beta/abs(phi1)*(2*eta/phi-meanL1)*(r-(1-exp(-r*abs(phi1)))/abs(phi1))+2*beta^2/phi^2*r)
#   phi2 <-2*termA*abs(phi1)/((var_G2-2*mu_G2^2)*abs(phi1)+termA)
#   if(typeacf == "covariance"){
#   TheoCovQuad <- meanL1*beta^2/abs(phi1)^3*(2*eta/phi-meanL1)*
#      (2/abs(phi2)-1/abs(phi1))*(1-exp(-r*abs(phi1)))*(exp(r*abs(phi1))-1)*
#      exp(-h*abs(phi1))
#   }else{
#     TheoCovQuad <- meanL1*beta^2/abs(phi1)^3*(2*eta/phi-meanL1)*
#       (2/abs(phi2)-1/abs(phi1))*(1-exp(-r*abs(phi1)))*(exp(r*abs(phi1))-1)*
#       exp(-h*abs(phi1))/(var_G2)
#   }
#   
#  } 
 if(env$q >= 1){
#    M2Lev <- param[env$M2Lev.idx] 
#    M4Lev <- param[env$M4Lev.idx]
#  }else{
#    M2Lev <- param[env$phi1.idx]
#    M4Lev <- param[env$phi2.idx]
#  }
   TheoCovQuad <- numeric(length = length(h))
   for(i in c(1:length(h))){
      MomentCog <- MM_Cogarch(p = env$p, q = env$q, cost=param[cost] , acoeff=param[a],
                              b=param[b], meanL1 = meanL1, r = r, h = h[i], 
                              type = typeacf, m2=mu_G2, var=var_G2)
   
   
   TheoCovQuad[i] <- MomentCog$acfG2
   }
   theo_mu_G2 <- MomentCog$meanG2
 }
 res <- sum((c(log(abs(TheoCovQuad)))-c(CovQuad))^2)
 return(res)
#  if(env$q == 3){
#    MomentCog <- MM_Cogarch(p = env$p, q = env$q, a, b, r, h)
#  }
#  if(env$q >= 4){
#    MomentCog <- MM_Cogarch(p = env$p, q = env$q, a, b, r, h)
#  }
}

MM_Cogarch <- function(p, q, cost, acoeff, b, meanL1, r, h, type, m2, var){
  # The code developed here is based on the acf for squared returns derived in the 
  # Chaadra phd Thesis
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  a[1:p,1] <- acoeff
#   mu <- M2Lev
  bq <- b[1]
  a1 <- a[1]
  mu <-1/acoeff[1]*(bq-cost*bq*r*meanL1/m2)
#   rho <- M4Lev
  B_tilde <- MatrixA(b[c(q:1)])+mu*e%*%t(a)
  #beta <- as.numeric(B_tilde[q,])
  meanG2 <- cost*bq*r/(bq-mu*a1)*meanL1
  Inf_eps <- IdentM <- diag(q)
  if(q==1){
    Inf_eps[q,q] <- 1/(-2*B_tilde[1,1]) 
  }
  if(q==2){
    Inf_eps[q,q] <- -B_tilde[2,1]
    Inf_eps <- 1/(2*B_tilde[2,1]*B_tilde[2,2])*Inf_eps
  }
  if(q==3){
    Inf_eps[1,1] <- B_tilde[3,3]/B_tilde[3,1]
    Inf_eps[q,q] <- -B_tilde[3,2]
    Inf_eps[1,3] <- -1
    Inf_eps[3,1] <- -1
    Inf_eps <- 1/(2*(B_tilde[3,3]*B_tilde[3,2]+B_tilde[3,1]))*Inf_eps
  }
  if(q>=4){
    Inf_eps <- round(AsympVar(B_tilde,e,lower=0,upper=100,delta=0.01)*10^5)/10^5
  }
   
  term <- expm(B_tilde*h)
  invB <- solve(B_tilde)
  term1 <- invB%*%(IdentM-expm(-B_tilde*r))
  term2 <- (expm(B_tilde*r)-IdentM)
  
  P0_overRho <- 2*mu^2*(3*invB%*%(invB%*%term2-r*IdentM)-IdentM)%*%Inf_eps
  Q0_overRho <- 6*mu*((r*IdentM-invB%*%term2)%*%Inf_eps
                      -invB%*%(invB%*%term2-r*IdentM)%*%Inf_eps%*%t(B_tilde))%*%e
  m_overRho <- as.numeric(t(a)%*%Inf_eps%*%a)
  Den<- (m_overRho*meanL1^2/m2^2*var*r^2+t(a)%*%Q0_overRho+t(a)%*%P0_overRho%*%a+1) 
  num <-(meanL1^2/m2^2*var-2*mu^2)*r^2
  rh0 <- as.numeric(num/Den)


  Inf_eps1 <- Inf_eps*rh0
  Ph <- mu^2*term%*%term1%*%invB%*%term2%*%Inf_eps1
  Qh <- mu*term%*%term1%*%(-term2%*%Inf_eps1-invB%*%term2%*%Inf_eps1%*%t(B_tilde))%*%e
  m <- m_overRho*rh0
  
  
 
  if(type=="correlation"){
    VarTheo<-as.numeric(rh0*t(a)%*%P0_overRho%*%a+rh0*t(a)%*%Q0_overRho+2*r^2*mu^2+rh0)
    acfG2 <- (t(a)%*%Ph%*%a+t(a)%*%Qh)/VarTheo
  }else{
    coeff <- cost^2*bq^2/((1-m)*(bq-mu*a1)^2)
    acfG2 <- coeff*(t(a)%*%Ph%*%a+t(a)%*%Qh)
  }
  res <- list(acfG2=acfG2, meanG2=meanG2)
  return(res)
}
AsympVar<-function(B_tilde,e,lower=0,upper=100,delta=0.1){
  part<-seq(lower,upper-delta, by=delta)+delta/2
  last <- length(part)
  Integrand <- as.matrix(matrix(0,length(e),length(e)))
  for(i in c(1:last)){
  Integrand<-Integrand+expm(B_tilde*part[i])%*%e%*%t(e)%*%expm(t(B_tilde)*part[i])*delta
  }
  return(Integrand)
}
