# In this code we implement a filter that returns the increments of the undelying lévy process
# if the model is a COGARCH(P,Q)
# Using the squared of returns, we obtain the increment of Levy process using the relation
# Y_t=e^{A\Delta t}Y_{t-\Delta t}+e^{A\left(\Delta t\right)}e\left(\Delta G_{t}\right)^{2}.
# for the state process $Y_t$ and then the levy increments definely are:
# \Delta L_{t}=\frac{\Delta G_{t}}{\sqrt{V_{t}}}.

cogarchNoise<-function(yuima.cogarch, data=NULL, param, mu=1){
  if(missing(yuima.cogarch))
    yuima.stop("yuima.cogarch or yuima object is missing.")
  
  if(length(param)==0)
    yuima.stop("missing values parameters")
  
  if(!is.COGARCH(yuima.cogarch)){
    yuima.warn("The model does not belong to the class yuima.cogarch")
  }
  
  if(is(yuima.cogarch,"yuima")){
    model<-yuima.cogarch@model
    if(is.null(data)){
      observ<-yuima.cogarch@data
    }else{observ<-data}
  }else{
    if(is(yuima.cogarch,"yuima.cogarch")){
      model<-yuima.cogarch
      if(is.null(data)){
        yuima.stop("Missing data")
      }
      observ<-data
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
  Data<-as.matrix(onezoo(observ)[,1])
  freq<-round(frequency(onezoo(observ)[,1]))
  res<-auxcogarch.noise(cost,b,acoeff,mu,Data,freq)
  return(res)
}
auxcogarch.noise<-function(cost,b,acoeff,mu,Data,freq){
  
  res<-StationaryMoments(cost,b,acoeff,mu)
  ExpY0<-res$ExpVar
  
  q<-length(b) 
  a <- e <- matrix(0,nrow=q,ncol=1)
  e[q,1] <- 1
  p<-length(acoeff)
  a[1:p,1] <- acoeff
  B <- MatrixA(b[c(q:1)])
  DeltaG <- diff(Data)
  squaredG <- DeltaG^2
  
  Process_Y <- res$ExpVar
  var_V<-cost + sum(acoeff*Process_Y)
  delta <- 1/freq
  for(t in c(2:(length(Data)))){  
    # Y_t=e^{A\Delta t}Y_{t-\Delta t}+e^{A\left(\Delta t\right)}e\left(\Delta G_{t}\right)^{2}
    Process_Y <- cbind(Process_Y, (expm(B*delta)%*%(Process_Y[,t-1]+e*squaredG[t-1])))
    var_V[t] <- cost + sum(acoeff*Process_Y[,t])
  }
  #\Delta L_{t}=\frac{\Delta G_{t}}{\sqrt{V_{t}}}.
  
  incr.L<- DeltaG/sqrt(var_V[c(2:(length(Data)))])
  return(incr.L)
}