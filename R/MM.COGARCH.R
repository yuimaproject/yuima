# In this function we consider different kind estimation procedures for COGARCH(P,Q)
# Model. 
is.COGARCH <- function(obj){
  if(is(obj,"yuima"))
    return(is(obj@model, "yuima.cogarch"))
  if(is(obj,"yuima.cogarch"))
    return(is(obj, "yuima.cogarch"))
  return(FALSE)
}
yuima.acf<-function(data,lag.max){
  mu<-mean(data)
  leng <-length(data)
  var1<-sum((data-mu)*(data-mu))/leng
  res1<-numeric(length=lag.max)
  for (t in 0:lag.max){
    h<-leng-lag.max
    res1[t+1]<-sum((data[(1+t):leng]-mu)*(data[1:h]-mu)) 
  }
  res1<-res1/(leng*var1) 
  acfr<-res1[2:(lag.max+1)]  #analogously to Matlab program
  
  
  
}

# The estimation procedure for cogarch(p,q) implemented in this code are based on the 
# Chadraa phd's thesis
gmm<-function(yuima, data = NULL, start, method="BFGS", fixed = list(), 
                           lower, upper, lag.max = NULL, aggr.G =TRUE){
  print <- FALSE

  call <- match.call()
  
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
  
  if( !missing(upper) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The upper requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }

  if( !missing(lower) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The lower constraints requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }

  if( !missing(fixed) && (method!="L-BFGS-B"||method!="brent")){
    yuima.warn("The fixed constraints requires L-BFGS-B or brent methods. We change method in  L-BFGS-B")
    method <- "L-BFGS-B"
  }
  
  # We identify the model parameters
  info <- model@info
  
  numb.ar <- info@q
  ar.name <- paste(info@ar.par,c(numb.ar:1),sep="")
  
  numb.ma <- info@p
  ma.name <- paste(info@ma.par,c(1:numb.ma),sep="")
  
  loc.par <- info@loc.par
  #xinit.name <- paste(info@Latent.var, c(1:numb.ar), sep = "")

  xinit.name0 <- model@parameter@xinit
  idx <- na.omit(match(c(loc.par, ma.name), xinit.name0))
  xinit.name <- xinit.name0[-idx]

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

  env <- new.env()
  n <- length(observ)[1]
  
  #Lag
  assign("lag", lag.max, envir=env)
  
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
  
  if(aggr.G==TRUE){
    #dt<-round(deltat(onezoo(observ)[,1])*10^5)/10^5
    NumbForUnit<-round(env$deltaData)
    G_i <- diff(env$Data[seq(1,length(env$Data),by=NumbForUnit)])
    r<-1
  }else{
    G_i <- diff(env$Data)  
    r <- 1/env$deltaData 
  }
  assign("G_i", G_i, envir=env)
  assign("r", r, envir=env)
  mu_G2 <- mean(G_i^2)
  assign("mu_G2", mu_G2, envir=env)
  var_G2 <- mean(G_i^4) - mu_G2^2
  assign("var_G2", var_G2, envir=env)
  d <- floor(sqrt(length(G_i)))
  assign("d", d, envir=env)
  typeacf <- "correlation"
  assign("typeacf", typeacf, envir=env)
  CovQuad <- log(abs(acf(G_i^2,plot=FALSE,lag.max=min(d,env$lag),type=typeacf)$acf[-1]))
  #CovQuad <-log(abs(yuima.acf(data=G_i^2,lag.max=min(d,env$lag))))
  assign("CovQuad", CovQuad, envir=env)

objectiveFun <- function(p,env) {
  mycoef <- as.list(p)
    if(length(c(fixed.idx, meas.idx))>0){ ## SMI 2/9/14
      names(mycoef) <- nm[-c(fixed.idx,meas.idx)] ## SMI 2/9/14
    }else{
      names(mycoef) <- nm
   } 
  ErrTerm(yuima=yuima, param=mycoef, print=print, env)
}
if(method!="L-BFGS-B"||method!="brent"){
  out<- optim(start, objectiveFun, method = method, env=env)
}else{
  if(length(fixed)!=0 && !missing(upper) && !missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                fixed=as.numeric(fixed), 
                lower=as.numeric(lower), 
                upper=as.numeric(upper), env=env)
  }else{
    if(!missing(upper) && !missing(lower)){
      out<- optim(start, objectiveFun, method = method, 
                  lower=as.numeric(lower), 
                  upper=as.numeric(upper), env=env)
    }
    if(length(fixed)!=0 && !missing(lower)){
      out<- optim(start, objectiveFun, method = method, 
                  fixed=as.numeric(fixed), 
                  lower=as.numeric(lower), env=env)
    }
    if(!missing(upper) && length(fixed)!=0){
      out<- optim(start, objectiveFun, method = method, 
                  fixed=as.numeric(fixed), 
                  upper=as.numeric(upper), env=env)
    }
  }
  
  if(length(fixed)!=0 && missing(upper) && missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                fixed=as.numeric(fixed), env=env)
  }
  
  if(length(fixed)==0 && !missing(upper) && missing(lower)){
    out<- optim(start, objectiveFun, method = method,
                lower=as.numeric(lower), env=env)
  }
  
  if(length(fixed)==0 && missing(upper) && !missing(lower)){
    out<- optim(start, objectiveFun, method = method, 
                upper=as.numeric(upper), env=env)
  }
  
}

 bvect<-out$par[ar.name]
 bq<-bvect[1]
 avect<-out$par[ma.name]
 a1<-avect[1]
 
  out$par[loc.par]<-(bq-a1)*mu_G2/(bq*r)

 # Determine the Variance-Covariance Matrix
 dimOutp<-length(out$par)-2
 coef <- out$par[c(1:dimOutp)] 
 vcov<-matrix(NA, dimOutp, dimOutp)
  names_coef<-names(coef)
  colnames(vcov) <- names_coef
  rownames(vcov) <- names_coef
  mycoef <- start
  min <- out$value 
# # call

 logL.Incr<-cogarchNoise(yuima.cogarch=model, data=observ, 
                    param=as.list(coef), mu=1)

# Build an object of class mle

 
 

 # Build an object of class cogarch.gmm.incr
res<-new("cogarch.gmm.incr", call = call, coef = coef, fullcoef = unlist(coef), 
    vcov = vcov, min = min, details = list(), 
    method = character(),
    Incr.Lev=logL.Incr,
    model = model, nobs=as.integer(length(logL.Incr)+1),
    logL.Incr = numeric())



 return(res)
  
}

ErrTerm <- function(yuima, param, print, env){
  typeacf <- env$typeacf
  param <- as.numeric(param)
  
  G_i <- env$G_i
  r <- env$r
  mu_G2 <- env$mu_G2
  var_G2 <- env$var_G2
  d <- env$d
  CovQuad <-env$CovQuad

  h <- seq(1, d, by = 1)*r
  cost <- env$loc.idx
  b <- env$ar.idx
  a <- env$ma.idx
  meanL1 <- param[env$EL1.idx]
#   meanL1 <- 1
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
# } 
 if(env$q >= 1){
   TheoCovQuad <- numeric(length = length(h))
   cost<-param[cost]
   for(i in c(1:length(h))){
      MomentCog <- MM_Cogarch(p = env$p, q = env$q,  acoeff=param[a],
                              b=param[b],  r = r, h = h[i], 
                              type = typeacf, m2=mu_G2, var=var_G2)
   
   
   TheoCovQuad[i] <- MomentCog$acfG2
   }
   theo_mu_G2 <- MomentCog$meanG2
   param[cost]<-MomentCog$cost
 }
 res <- sum((c(log(abs(TheoCovQuad)))-c(CovQuad))^2)
 return(res)
}

MM_Cogarch <- function(p, q, acoeff, b,  r, h, type, m2, var){
  # The code developed here is based on the acf for squared returns derived in the 
  # Chaadra phd Thesis
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  a[1:p,1] <- acoeff

  bq <- b[1]
  a1 <- a[1]

# # Matching only the autocorrelation we are not able to estimate the a_0 parameter
# nevertheless using the theoretical formula of the expectaion of squared returns we
# are able to fix this parameter for having a perfect match  between theoretical and 
# empirical mean

# We recall that under the assumption of levy  process is centered  the mean at time 1 of the 
# squared returns and the mean of the corresponding levy measure are equals.

  mu<-1 # we assume variance of the underlying L\'evy is one
  meanL1<-mu

  cost<-(bq-mu*a1)*m2/(bq*r)
  
  B_tilde <- MatrixA(b[c(q:1)])+mu*e%*%t(a)
  meanG2 <- cost*bq*r/(bq-mu*a1)*meanL1
  Inf_eps <- IdentM <- diag(q)
  if(q==1){
    Inf_eps[q,q] <- -1/(2*B_tilde[1,1]) 
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
  invB <- solve(B_tilde) # In this case we can use the analytical form for Companion Matrix???
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
  res <- list(acfG2=acfG2, meanG2=meanG2, cost=cost)
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
