yuima.PhamBreton.Alg<-function(a){
  p<-length(a)
  gamma<-a[p:1]
  if(p>2){
    gamma[p]<-a[1]
    alpha<-matrix(NA,p,p)
    for(j in 1:p){
      if(is.integer(as.integer(j)/2)){
        alpha[p,j]<-0
        alpha[p-1,j]<-0
      }else{
        alpha[p,j]<-a[j]
        alpha[p-1,j]<-a[j+1]/gamma[p]
      }
    }
    for(n in (p-1):1){
      gamma[n]<-alpha[n+1,2]-alpha[n,2]
      for(j in 1:n-1){
        alpha[n-1,j]<-(alpha[n+1,j+2]-alpha[n,j+2])/gamma[n]
      }
      alpha[n-1,n-1]<-alpha[n+1,n+1]/gamma[n]
    }
    gamma[1]<-alpha[2,2]
  }
  return(gamma)
}


Diagnostic.Carma<-function(carma){

  if(!is(carma@model,"yuima.carma"))
    yuima.stop("model is not a carma")
  if(!is(carma,"mle"))
    yuima.stop("object does not belong
                to yuima.qmle-class or yuima.carma.qmle-class ")
  param<-coef(carma)
  info <- carma@model@info
  numb.ar<-info@p
  name.ar<-paste(info@ar.par,c(numb.ar:1),sep="")
  ar.par<-param[name.ar]

  statCond<-FALSE
  if(min(yuima.PhamBreton.Alg(ar.par[numb.ar:1]))>=0)
    statCond<-TRUE
  return(statCond)
}

