
# In this file we develop the procedure described in Brockwell, Davis and Yang (2012) 
# for estimation of the underlying noise once the parameters of carma(p,q) are previously 
# obtained
 

yuima.carma.eigen<-function(A){
  diagA<-eigen(A)
  diagA$values<-diagA$values[order(diagA$values, na.last = TRUE, decreasing = TRUE)]
  n_eigenval<-length(diagA$values)
  diagA$vectors<-matrix(diagA$values[1]^(c(1:n_eigenval)-1),n_eigenval,1)
  if(n_eigenval>=2){
    for (i in 2:n_eigenval){
      diagA$vectors<-cbind(diagA$vectors, 
                           matrix(diagA$values[i]^(c(1:n_eigenval)-1),n_eigenval,1))
    }
  }
  return(diagA)
}


StateVarX<-function(y,tt,X_0,B){
  # The code obtains the first q-1 state variable using eq 5.1 in Brockwell, Davis and Yang 2011
  Time<-length(tt) 
  q<-length(X_0)
  e_q<-rep(0,q)
  e_q[q]<-1
  X_0<-matrix(X_0,q,1)
  e_q<-matrix(e_q,q,1)
  X<-matrix(0,q,Time)
  int<-matrix(0,q,1)
  for (i in c(2:Time)){
    int<-int+(expm(B*(tt[i]-tt[(i-1)]))%*%(e_q*y[i]))*(tt[i]-tt[(i-1)])
    X[i]<-as.matrix(expm(B*tt[i])%*%X_0+int)
  }
  return(X)
}


StateVarXp<-function(y,X_q,tt,B,q,p){
  # The code computes the  state variable X using the derivatives of eq 5.2 
  # see Brockwell, Davis and Yang 2011
  
  diagMatB<-yuima.carma.eigen(B)
  if(length(diagMatB$values)>1){
    MatrD<-diag(diagMatB$values)
  }else{
    MatrD<-as.matrix(diagMatB$values)
  }
  MatrR<-diagMatB$vectors
  idx.r<-c(q:(p-1))
  elem.X <-length(idx.r)
  YMatr<-matrix(0,q,length(y))
  YMatr[q,]<-y
  OutherX<-matrix(NA,q,length(y))
  OutherXalt<-matrix(NA,q,length(y))
  for(i in 1:elem.X){
    OutherX[i,]<-((MatrR%*%MatrD^(idx.r[i])%*%solve(MatrR))%*%X_q+
      (MatrR%*%MatrD^(idx.r[i]-1)%*%solve(MatrR))%*%YMatr)[1,]  
  }
  X.StatVar<-rbind(X_q,OutherX)
  return(X.StatVar)
}

bEvalPoly<-function(b,lambdax){
  result<-sum(b*lambdax^(c(1:length(b))-1))
  return(result)
}

aEvalPoly<-function(a,lambdax){
  p<-length(a)
  a.new<-c(1,a[c(1:(p-1))])
  pa.new<-c(p:1)*a.new
  result<-sum(pa.new*lambdax^(c(p:1)-1))
  return(result)
}


CarmaRecovNoise<-function(yuima, param, data=NULL){
  if( missing(param) ) 
    yuima.stop("Parameter values are missing.")
  
  if(!is.list(param))
    yuima.stop("Argument 'param' must be of list type.")
  
  vect.param<-as.numeric(param)
  name.param<-names(param)
  names(vect.param)<-name.param
  
  if(is(yuima,"yuima")){
    model<-yuima@model
    if(is.null(data)){
      observ<-yuima@data
    }else{observ<-data}
}else{
  if(is(yuima,"yuima.carma")){
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
  
  info<-model@info
  
  numb.ar<-info@p
  name.ar<-paste(info@ar.par,c(numb.ar:1),sep="")
  ar.par<-vect.param[name.ar]
  
  numb.ma<-info@q
  name.ma<-paste(info@ma.par,c(0:numb.ma),sep="")
  ma.par<-vect.param[name.ma]
  
  loc.par=NULL
  if (length(info@loc.par)!=0){
    loc.par<-vect.param[info@loc.par]
  }
  
  scale.par=NULL
  if (length(info@scale.par)!=0){
    scale.par<-vect.param[info@scale.par]
  }
  
  lin.par=NULL
  if (length(info@lin.par)!=0){
    lin.par<-vect.param[info@lin.par]
  }
  
  
  
  ttt<-observ@zoo.data[[1]]
  tt<-index(ttt)
  y<-coredata(ttt)
  
  levy<-yuima.CarmaRecovNoise(y,tt,ar.par,ma.par, loc.par, scale.par, lin.par)
  inc.levy<-diff(t(levy))
  return(inc.levy)
}



yuima.CarmaRecovNoise<-function(y,tt,ar.par,ma.par, loc.par=NULL, scale.par=NULL, lin.par=NULL){
    if(!is.null(loc.par)){
      yuima.warn("the loc.par will be implemented as soon as possible")
      return(NULL)
    } 
    if(!is.null(scale.par)){
      yuima.warn("the scale.par will be implemented as soon as possible")
      return(NULL)
    }
    if(!is.null(lin.par)){
      yuima.warn("the lin.par will be implemented as soon as possible")
      return(NULL)
    }
    
    p<-length(ar.par)
    q<-length(ma.par)
    
    if(q==1){
      yuima.warn("the car(p) process will be implemented as soon as possible")
      return(NULL)
    }else{
    A<-MatrixA(ar.par[c(p:1)])
    
    b_q<-tail(ma.par,n=1)
    
    newma.par<-ma.par/b_q
    # We build the matrix B which is necessary for building eq 5.2
    ynew<-y/b_q
    B<-MatrixA(newma.par[c(1:(q-1))])
    diagB<-yuima.carma.eigen(B)
    
    e_q<-rep(0,(q-1))
    if((q-1)>0){
      e_q[(q-1)]<-1
      X_0<-rep(0,(q-1))
    }else{
      e_q<-1
      X_0<-0
    }
    
    
 
    X_q<-StateVarX(ynew,tt,X_0,B)
    
    #plot(t(X_q))
    
    X.StVa<-StateVarXp(ynew,X_q,tt,B,q-1,p)
    
    #plot(y)
    
    diagA<-yuima.carma.eigen(A)

  
       
    
    BinLambda<-rep(NA,length(ma.par))
    for(i in c(1:length(diagA$values))){
      BinLambda[i]<-bEvalPoly(ma.par,diagA$values[i])
    }
    MatrBLam<-diag(BinLambda)
    
    # We get the Canonical Vector Space in eq 2.17
    Y_CVS<-MatrBLam%*%solve(diagA$vectors)%*%X.StVa #Canonical Vector Space
  
    # We verify the prop 2 in the paper "Estimation for Non-Negative Levy Driven CARMA process
    # yver<-Y_CVS[1,]+Y_CVS[2,]
    
    # plot(yver)
    
    # plot(y)
    
    
    idx.r<-match(0,Im(diagA$values))
    lambda.r<-diagA$values[idx.r]
    int<-0
    lev.und<-matrix(0,1,length(tt))
    derA<-aEvalPoly(ar.par[c(p:1)],lambda.r)
      for(t in c(2:length(tt))){
        int<-int+y[t]*(tt[t]-tt[t-1])
        lev.und[,t]<-derA/BinLambda[idx.r]*(Y_CVS[idx.r,t]-Y_CVS[idx.r,1]-lambda.r*int)
      }
    }
  return(lev.und)
}

